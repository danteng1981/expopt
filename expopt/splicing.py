"""
Splicing Module

This module provides functionality for:
- RGroup Splicing: Connecting replacement fragments at attachment points via single bonds
- Core Splicing: Enumerating connection point combinations and scoring by fingerprint similarity
"""

from typing import List, Optional, Dict, Any, Tuple, Iterator
from dataclasses import dataclass
from itertools import permutations, product

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
    from rdkit.Chem.Fingerprints import FingerprintMols
    from rdkit import DataStructs
    from rdkit.Chem import rdmolops
except ImportError:
    Chem = None
    AllChem = None
    rdMolDescriptors = None
    FingerprintMols = None
    DataStructs = None
    rdmolops = None


@dataclass
class SplicedMolecule:
    """Represents a molecule created by splicing."""
    smiles: str
    mol: Any  # RDKit Mol object
    similarity_to_target: float = 0.0
    connection_info: Optional[Dict[str, Any]] = None

    def __post_init__(self):
        if self.connection_info is None:
            self.connection_info = {}


class RGroupSplicer:
    """
    RGroup Splicing Engine

    Connects replacement fragments at specified attachment points
    via single bonds to the original scaffold.
    """

    def __init__(self):
        if Chem is None:
            raise ImportError("RDKit is required for this functionality")

    def splice(
        self,
        scaffold_smiles: str,
        fragment_smiles: str,
        scaffold_attach_idx: int,
        fragment_attach_idx: int
    ) -> Optional[SplicedMolecule]:
        """
        Splice a fragment onto a scaffold at specified attachment points.

        Args:
            scaffold_smiles: SMILES of the scaffold molecule
            fragment_smiles: SMILES of the fragment to attach
            scaffold_attach_idx: Atom index on scaffold for connection
            fragment_attach_idx: Atom index on fragment for connection

        Returns:
            SplicedMolecule if successful, None otherwise
        """
        scaffold_mol = Chem.MolFromSmiles(scaffold_smiles)
        fragment_mol = Chem.MolFromSmiles(fragment_smiles)

        if scaffold_mol is None or fragment_mol is None:
            return None

        try:
            # Combine molecules
            combined = Chem.CombineMols(scaffold_mol, fragment_mol)
            editable = Chem.EditableMol(combined)

            # Adjust fragment atom index (offset by scaffold size)
            offset_fragment_idx = fragment_attach_idx + scaffold_mol.GetNumAtoms()

            # Add single bond between attachment points
            editable.AddBond(
                scaffold_attach_idx,
                offset_fragment_idx,
                Chem.BondType.SINGLE
            )

            # Get the result
            result_mol = editable.GetMol()
            Chem.SanitizeMol(result_mol)

            return SplicedMolecule(
                smiles=Chem.MolToSmiles(result_mol),
                mol=result_mol,
                connection_info={
                    "scaffold_idx": scaffold_attach_idx,
                    "fragment_idx": fragment_attach_idx,
                }
            )
        except Exception:
            return None

    def splice_with_dummy_atoms(
        self,
        scaffold_smiles: str,
        fragment_smiles: str
    ) -> Optional[SplicedMolecule]:
        """
        Splice molecules using dummy atoms (*) as attachment points.

        The dummy atoms are removed and replaced with a single bond
        connecting the real atoms.

        Args:
            scaffold_smiles: SMILES with dummy atom marking attachment
            fragment_smiles: SMILES with dummy atom marking attachment

        Returns:
            SplicedMolecule if successful, None otherwise
        """
        scaffold_mol = Chem.MolFromSmiles(scaffold_smiles)
        fragment_mol = Chem.MolFromSmiles(fragment_smiles)

        if scaffold_mol is None or fragment_mol is None:
            return None

        # Find dummy atoms and their neighbors
        scaffold_attach = self._find_dummy_attachment(scaffold_mol)
        fragment_attach = self._find_dummy_attachment(fragment_mol)

        if scaffold_attach is None or fragment_attach is None:
            return None

        scaffold_dummy_idx, scaffold_real_idx = scaffold_attach
        fragment_dummy_idx, fragment_real_idx = fragment_attach

        # Remove dummy atoms and connect real atoms
        return self._connect_after_removing_dummies(
            scaffold_mol, fragment_mol,
            scaffold_dummy_idx, scaffold_real_idx,
            fragment_dummy_idx, fragment_real_idx
        )

    def _find_dummy_attachment(self, mol: Any) -> Optional[Tuple[int, int]]:
        """Find dummy atom and its connected real atom."""
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:  # Dummy atom
                neighbors = atom.GetNeighbors()
                if len(neighbors) == 1:
                    return (atom.GetIdx(), neighbors[0].GetIdx())
        return None

    def _connect_after_removing_dummies(
        self,
        scaffold_mol: Any,
        fragment_mol: Any,
        scaffold_dummy_idx: int,
        scaffold_real_idx: int,
        fragment_dummy_idx: int,
        fragment_real_idx: int
    ) -> Optional[SplicedMolecule]:
        """Connect molecules after removing dummy atoms."""
        try:
            # Remove dummy atoms from copies
            scaffold_edit = Chem.RWMol(scaffold_mol)
            scaffold_edit.RemoveAtom(scaffold_dummy_idx)

            fragment_edit = Chem.RWMol(fragment_mol)
            fragment_edit.RemoveAtom(fragment_dummy_idx)

            # Adjust indices after removal
            new_scaffold_idx = scaffold_real_idx
            if scaffold_real_idx > scaffold_dummy_idx:
                new_scaffold_idx -= 1

            new_fragment_idx = fragment_real_idx
            if fragment_real_idx > fragment_dummy_idx:
                new_fragment_idx -= 1

            # Combine and connect
            combined = Chem.CombineMols(scaffold_edit, fragment_edit)
            editable = Chem.EditableMol(combined)

            offset_fragment_idx = new_fragment_idx + scaffold_edit.GetNumAtoms()
            editable.AddBond(
                new_scaffold_idx,
                offset_fragment_idx,
                Chem.BondType.SINGLE
            )

            result_mol = editable.GetMol()
            Chem.SanitizeMol(result_mol)

            return SplicedMolecule(
                smiles=Chem.MolToSmiles(result_mol),
                mol=result_mol,
                connection_info={
                    "scaffold_idx": new_scaffold_idx,
                    "fragment_idx": new_fragment_idx,
                }
            )
        except Exception:
            return None

    def splice_multiple_fragments(
        self,
        scaffold_smiles: str,
        fragments: List[Tuple[str, int, int]]
    ) -> Optional[SplicedMolecule]:
        """
        Splice multiple fragments onto a scaffold.

        Args:
            scaffold_smiles: SMILES of the scaffold
            fragments: List of (fragment_smiles, scaffold_idx, fragment_idx) tuples

        Returns:
            SplicedMolecule with all fragments attached
        """
        current_mol = Chem.MolFromSmiles(scaffold_smiles)
        if current_mol is None:
            return None

        connection_info = {"fragments": []}
        atom_offset = 0

        for frag_smiles, scaffold_idx, frag_idx in fragments:
            frag_mol = Chem.MolFromSmiles(frag_smiles)
            if frag_mol is None:
                continue

            try:
                combined = Chem.CombineMols(current_mol, frag_mol)
                editable = Chem.EditableMol(combined)

                offset_frag_idx = frag_idx + current_mol.GetNumAtoms()
                adjusted_scaffold_idx = scaffold_idx + atom_offset

                editable.AddBond(
                    adjusted_scaffold_idx,
                    offset_frag_idx,
                    Chem.BondType.SINGLE
                )

                current_mol = editable.GetMol()
                Chem.SanitizeMol(current_mol)

                connection_info["fragments"].append({
                    "smiles": frag_smiles,
                    "scaffold_idx": adjusted_scaffold_idx,
                    "fragment_idx": frag_idx,
                })
            except Exception:
                continue

        return SplicedMolecule(
            smiles=Chem.MolToSmiles(current_mol),
            mol=current_mol,
            connection_info=connection_info
        )


class CoreSplicer:
    """
    Core Splicing Engine

    Enumerates all reasonable connection point combinations between
    candidate cores and side chains/R-groups, scores by fingerprint
    similarity, and selects the most promising connections.
    """

    def __init__(self):
        if Chem is None:
            raise ImportError("RDKit is required for this functionality")

    def enumerate_connections(
        self,
        core_smiles: str,
        rgroups: List[str],
        core_attachment_points: Optional[List[int]] = None
    ) -> Iterator[Tuple[List[Tuple[int, int]], str]]:
        """
        Enumerate all possible connection combinations between a core
        and R-groups.

        Args:
            core_smiles: SMILES of the core scaffold
            rgroups: List of R-group SMILES
            core_attachment_points: Specific attachment points on core (optional)

        Yields:
            Tuples of (connection_list, description) where connection_list
            contains (core_idx, rgroup_idx) pairs
        """
        core_mol = Chem.MolFromSmiles(core_smiles)
        if core_mol is None:
            return

        # Get attachment points if not specified
        if core_attachment_points is None:
            core_attachment_points = self._find_attachment_points(core_mol)

        if len(core_attachment_points) < len(rgroups):
            return

        # Get attachment points for each R-group
        rgroup_attachments = []
        for rg_smiles in rgroups:
            rg_mol = Chem.MolFromSmiles(rg_smiles)
            if rg_mol is None:
                return
            rg_attach = self._find_attachment_points(rg_mol)
            if not rg_attach:
                rg_attach = [0]  # Default to first atom
            rgroup_attachments.append(rg_attach)

        # Enumerate permutations of core attachment points
        for core_perm in permutations(core_attachment_points, len(rgroups)):
            # For each permutation, enumerate R-group attachment combinations
            for rg_combo in product(*rgroup_attachments):
                connections = list(zip(core_perm, rg_combo))
                desc = f"Core attachments: {core_perm}, RGroup attachments: {rg_combo}"
                yield connections, desc

    def splice_core_with_rgroups(
        self,
        core_smiles: str,
        rgroups: List[str],
        connections: List[Tuple[int, int]]
    ) -> Optional[SplicedMolecule]:
        """
        Splice a core with R-groups using specified connections.

        Args:
            core_smiles: SMILES of the core scaffold
            rgroups: List of R-group SMILES
            connections: List of (core_attach_idx, rgroup_attach_idx) tuples

        Returns:
            SplicedMolecule if successful, None otherwise
        """
        if len(rgroups) != len(connections):
            return None

        core_mol = Chem.MolFromSmiles(core_smiles)
        if core_mol is None:
            return None

        current_mol = core_mol
        atom_offset = 0
        connection_info = {"core": core_smiles, "rgroups": [], "bonds": []}

        for i, (rg_smiles, (core_idx, rg_idx)) in enumerate(zip(rgroups, connections)):
            rg_mol = Chem.MolFromSmiles(rg_smiles)
            if rg_mol is None:
                continue

            try:
                combined = Chem.CombineMols(current_mol, rg_mol)
                editable = Chem.EditableMol(combined)

                # Calculate actual indices
                actual_core_idx = core_idx  # Core indices don't change
                actual_rg_idx = rg_idx + current_mol.GetNumAtoms()

                editable.AddBond(
                    actual_core_idx,
                    actual_rg_idx,
                    Chem.BondType.SINGLE
                )

                current_mol = editable.GetMol()
                Chem.SanitizeMol(current_mol)

                connection_info["rgroups"].append(rg_smiles)
                connection_info["bonds"].append({
                    "core_idx": core_idx,
                    "rgroup_idx": rg_idx,
                })

                atom_offset += rg_mol.GetNumAtoms()
            except Exception:
                return None

        return SplicedMolecule(
            smiles=Chem.MolToSmiles(current_mol),
            mol=current_mol,
            connection_info=connection_info
        )

    def find_best_connections(
        self,
        core_smiles: str,
        rgroups: List[str],
        target_smiles: str,
        max_results: int = 10,
        core_attachment_points: Optional[List[int]] = None
    ) -> List[SplicedMolecule]:
        """
        Find the best connection combinations by scoring against target molecule.

        For each candidate core, enumerate all reasonable connection point
        combinations, create complete molecules, and score by fingerprint
        similarity to the target molecule.

        Args:
            core_smiles: SMILES of the core scaffold
            rgroups: List of R-group SMILES
            target_smiles: SMILES of the target molecule for similarity scoring
            max_results: Maximum number of results to return
            core_attachment_points: Specific attachment points on core

        Returns:
            List of SplicedMolecule objects sorted by similarity score
        """
        target_mol = Chem.MolFromSmiles(target_smiles)
        if target_mol is None:
            return []

        target_fp = FingerprintMols.FingerprintMol(target_mol)
        results = []

        for connections, _ in self.enumerate_connections(
            core_smiles, rgroups, core_attachment_points
        ):
            spliced = self.splice_core_with_rgroups(core_smiles, rgroups, connections)
            if spliced is None or spliced.mol is None:
                continue

            # Calculate similarity to target
            spliced_fp = FingerprintMols.FingerprintMol(spliced.mol)
            similarity = DataStructs.TanimotoSimilarity(target_fp, spliced_fp)
            spliced.similarity_to_target = similarity

            results.append(spliced)

        # Sort by similarity and return top results
        results.sort(key=lambda x: x.similarity_to_target, reverse=True)
        return results[:max_results]

    def _find_attachment_points(self, mol: Any) -> List[int]:
        """
        Find potential attachment points on a molecule.

        Looks for atoms with available valence (implicit hydrogens)
        or explicit dummy atoms.
        """
        attachment_points = []

        for atom in mol.GetAtoms():
            # Dummy atoms are explicit attachment points
            if atom.GetAtomicNum() == 0:
                attachment_points.append(atom.GetIdx())
                continue

            # Atoms with hydrogens could be attachment points
            if atom.GetTotalNumHs() > 0:
                attachment_points.append(atom.GetIdx())

        return attachment_points

    def score_molecule(
        self,
        mol_smiles: str,
        target_smiles: str
    ) -> float:
        """
        Calculate Tanimoto fingerprint similarity between two molecules.

        Args:
            mol_smiles: SMILES of the molecule to score
            target_smiles: SMILES of the target molecule

        Returns:
            Similarity score (0.0 to 1.0)
        """
        mol = Chem.MolFromSmiles(mol_smiles)
        target = Chem.MolFromSmiles(target_smiles)

        if mol is None or target is None:
            return 0.0

        mol_fp = FingerprintMols.FingerprintMol(mol)
        target_fp = FingerprintMols.FingerprintMol(target)

        return DataStructs.TanimotoSimilarity(mol_fp, target_fp)

    def batch_splice_and_score(
        self,
        cores: List[str],
        rgroups: List[str],
        target_smiles: str,
        top_n: int = 10
    ) -> List[SplicedMolecule]:
        """
        Process multiple candidate cores, splice with R-groups,
        and return the top scoring combinations.

        Args:
            cores: List of candidate core SMILES
            rgroups: List of R-group SMILES
            target_smiles: Target molecule for similarity scoring
            top_n: Number of top results to return

        Returns:
            List of top-scoring SplicedMolecule objects
        """
        all_results = []

        for core_smiles in cores:
            results = self.find_best_connections(
                core_smiles, rgroups, target_smiles,
                max_results=top_n
            )
            all_results.extend(results)

        # Sort all results and return top N
        all_results.sort(key=lambda x: x.similarity_to_target, reverse=True)
        return all_results[:top_n]
