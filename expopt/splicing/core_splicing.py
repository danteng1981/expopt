"""Core splicing implementation for scaffold replacement."""

from dataclasses import dataclass
from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem

from expopt.utils.similarity import SimilarityCalculator


@dataclass
class SplicedMolecule:
    """Represents a spliced molecule with similarity score."""

    smiles: str
    mol: Chem.Mol
    similarity_score: float
    configuration_id: int


class CoreSplicing:
    """
    Splice sidechains onto candidate cores.

    Enumerates connection configurations and ranks complete
    molecules by fingerprint similarity.
    """

    def __init__(self, fingerprint_type: str = "morgan"):
        """
        Initialize core splicing utility.

        Args:
            fingerprint_type: Type of fingerprint for similarity
        """
        self.similarity = SimilarityCalculator(fingerprint_type=fingerprint_type)

    def splice_core(
        self,
        core_smiles: str,
        sidechains: list[str],
        attachment_map: Optional[dict[int, int]] = None,
    ) -> Optional[str]:
        """
        Splice sidechains onto a core at specified positions.

        Args:
            core_smiles: SMILES of core scaffold
            sidechains: List of sidechain SMILES
            attachment_map: Dict mapping sidechain index to core atom index

        Returns:
            SMILES of spliced molecule or None
        """
        core_mol = Chem.MolFromSmiles(core_smiles)
        if core_mol is None:
            return None

        # If no attachment map, use default positions
        if attachment_map is None:
            attachment_map = {}
            available_positions = self._get_attachment_positions(core_mol)
            for i, _ in enumerate(sidechains):
                if i < len(available_positions):
                    attachment_map[i] = available_positions[i]

        # Build combined molecule
        combined = Chem.RWMol(core_mol)
        offset = core_mol.GetNumAtoms()

        sidechain_mols = []
        for sc_smiles in sidechains:
            sc_mol = Chem.MolFromSmiles(sc_smiles)
            if sc_mol is not None:
                sidechain_mols.append(sc_mol)

        # Add sidechains to combined molecule
        for i, sc_mol in enumerate(sidechain_mols):
            combined = Chem.RWMol(Chem.CombineMols(combined.GetMol(), sc_mol))

            if i in attachment_map:
                core_attach = attachment_map[i]
                # First atom of sidechain after offset
                sc_attach = offset

                try:
                    combined.AddBond(core_attach, sc_attach, Chem.BondType.SINGLE)
                except Exception:
                    pass

            offset += sc_mol.GetNumAtoms()

        try:
            result = combined.GetMol()
            Chem.SanitizeMol(result)
            return Chem.MolToSmiles(result)
        except Exception:
            return None

    def _get_attachment_positions(self, mol: Chem.Mol) -> list[int]:
        """
        Get potential attachment positions on a molecule.

        Args:
            mol: RDKit molecule

        Returns:
            List of atom indices suitable for attachment
        """
        positions = []

        for atom in mol.GetAtoms():
            # Check for dummy atoms first
            if atom.GetAtomicNum() == 0:
                neighbors = atom.GetNeighbors()
                if neighbors:
                    positions.append(neighbors[0].GetIdx())
            # Check for atoms with available valence
            elif self._has_available_valence(atom):
                positions.append(atom.GetIdx())

        return positions

    def _has_available_valence(self, atom: Chem.Atom) -> bool:
        """Check if atom has available valence for bonding."""
        total_valence = atom.GetTotalValence()
        default_valence = Chem.GetPeriodicTable().GetDefaultValence(atom.GetAtomicNum())

        if isinstance(default_valence, tuple):
            max_valence = max(default_valence)
        else:
            max_valence = default_valence

        return total_valence < max_valence

    def enumerate_configurations(
        self,
        core_smiles: str,
        sidechains: list[str],
        max_configurations: int = 100,
    ) -> list[SplicedMolecule]:
        """
        Enumerate different connection configurations.

        Args:
            core_smiles: SMILES of core
            sidechains: List of sidechain SMILES
            max_configurations: Maximum number of configurations to generate

        Returns:
            List of SplicedMolecule objects
        """
        core_mol = Chem.MolFromSmiles(core_smiles)
        if core_mol is None:
            return []

        positions = self._get_attachment_positions(core_mol)
        if len(positions) < len(sidechains):
            return []

        configurations = []
        config_id = 0
        seen_smiles = set()

        # Generate permutations of attachment positions
        from itertools import permutations

        for perm in permutations(positions, len(sidechains)):
            if config_id >= max_configurations:
                break

            attachment_map = {i: pos for i, pos in enumerate(perm)}
            result = self.splice_core(core_smiles, sidechains, attachment_map)

            if result and result not in seen_smiles:
                mol = Chem.MolFromSmiles(result)
                if mol is not None:
                    configurations.append(
                        SplicedMolecule(
                            smiles=result,
                            mol=mol,
                            similarity_score=0.0,  # Will be set during ranking
                            configuration_id=config_id,
                        )
                    )
                    seen_smiles.add(result)
                    config_id += 1

        return configurations

    def rank_by_similarity(
        self,
        candidates: list[SplicedMolecule],
        reference_mol: Chem.Mol,
    ) -> list[SplicedMolecule]:
        """
        Rank spliced molecules by fingerprint similarity.

        Args:
            candidates: List of spliced molecules
            reference_mol: Reference molecule for comparison

        Returns:
            Ranked list of candidates
        """
        for candidate in candidates:
            candidate.similarity_score = self.similarity.tanimoto_similarity(
                candidate.mol, reference_mol
            )

        # Sort by similarity score descending
        candidates.sort(key=lambda c: c.similarity_score, reverse=True)
        return candidates

    def splice_and_rank(
        self,
        core_smiles: str,
        sidechains: list[str],
        reference_smiles: str,
        max_configurations: int = 100,
        top_n: int = 10,
    ) -> list[SplicedMolecule]:
        """
        Complete workflow: enumerate configurations and rank by similarity.

        Args:
            core_smiles: Core scaffold SMILES
            sidechains: List of sidechain SMILES
            reference_smiles: Reference molecule SMILES for similarity
            max_configurations: Maximum configurations to generate
            top_n: Number of top results to return

        Returns:
            Top ranked spliced molecules
        """
        reference_mol = Chem.MolFromSmiles(reference_smiles)
        if reference_mol is None:
            return []

        candidates = self.enumerate_configurations(
            core_smiles, sidechains, max_configurations
        )

        if not candidates:
            return []

        ranked = self.rank_by_similarity(candidates, reference_mol)
        return ranked[:top_n]

    def get_best_configuration(
        self,
        core_smiles: str,
        sidechains: list[str],
        reference_smiles: str,
    ) -> Optional[str]:
        """
        Get the best configuration by similarity.

        Args:
            core_smiles: Core scaffold SMILES
            sidechains: Sidechain SMILES list
            reference_smiles: Reference molecule SMILES

        Returns:
            SMILES of best configuration or None
        """
        results = self.splice_and_rank(
            core_smiles, sidechains, reference_smiles, top_n=1
        )

        if results:
            return results[0].smiles
        return None
