"""
Core Hopping Module

This module provides functionality for replacing the core scaffold of molecules
while preserving side chains and key functional groups.

The goal is to find scaffolds with "hopping" potential - new core structures
that maintain the pharmacophoric features while providing novel chemical space.
"""

from typing import List, Optional, Dict, Any, Tuple, Set
from dataclasses import dataclass

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Chem.Fingerprints import FingerprintMols
    from rdkit import DataStructs
except ImportError:
    Chem = None
    AllChem = None
    Descriptors = None
    MurckoScaffold = None
    FingerprintMols = None
    DataStructs = None


@dataclass
class CoreCandidate:
    """Represents a candidate core scaffold for hopping."""
    smiles: str
    mol: Any  # RDKit Mol object
    similarity_score: float = 0.0
    attachment_points: List[int] = None
    properties: Optional[Dict[str, float]] = None
    source: str = "generated"

    def __post_init__(self):
        if self.attachment_points is None:
            self.attachment_points = []
        if self.properties is None:
            self.properties = {}


class CoreDatabase:
    """Database of core scaffolds for hopping."""

    def __init__(self):
        self._cores: Dict[str, Any] = {}  # canonical SMILES -> Mol object
        self._fingerprints: Dict[str, Any] = {}

    def add_core(self, smiles: str) -> bool:
        """Add a core scaffold to the database."""
        if Chem is None:
            raise ImportError("RDKit is required for this functionality")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        canonical_smiles = Chem.MolToSmiles(mol)
        self._cores[canonical_smiles] = mol
        self._fingerprints[canonical_smiles] = FingerprintMols.FingerprintMol(mol)
        return True

    def add_cores_from_list(self, smiles_list: List[str]) -> int:
        """Add multiple cores to the database. Returns count of successfully added."""
        count = 0
        for smiles in smiles_list:
            if self.add_core(smiles):
                count += 1
        return count

    def lookup(self, smiles: str) -> Optional[Any]:
        """Look up a core by SMILES. Returns Mol object if found."""
        if Chem is None:
            raise ImportError("RDKit is required for this functionality")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        canonical_smiles = Chem.MolToSmiles(mol)
        return self._cores.get(canonical_smiles)

    def find_similar(
        self, smiles: str, threshold: float = 0.3, top_n: int = 20
    ) -> List[Tuple[str, float]]:
        """Find cores similar to the query."""
        if Chem is None:
            raise ImportError("RDKit is required for this functionality")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []

        query_fp = FingerprintMols.FingerprintMol(mol)
        results = []

        for core_smiles, core_fp in self._fingerprints.items():
            similarity = DataStructs.TanimotoSimilarity(query_fp, core_fp)
            if similarity >= threshold:
                results.append((core_smiles, similarity))

        results.sort(key=lambda x: x[1], reverse=True)
        return results[:top_n]

    def get_all_cores(self) -> List[str]:
        """Return all core SMILES in the database."""
        return list(self._cores.keys())

    def __len__(self) -> int:
        return len(self._cores)


class CoreHopper:
    """
    Core Hopping Engine

    Replaces the core scaffold of molecules while preserving side chains
    and key functional groups to find scaffolds with hopping potential.
    """

    def __init__(self, database: Optional[CoreDatabase] = None):
        """
        Initialize the Core Hopper.

        Args:
            database: Optional core database for scaffold matching
        """
        if Chem is None:
            raise ImportError("RDKit is required for this functionality")

        self.database = database or CoreDatabase()

    def extract_core(self, mol_smiles: str) -> Optional[str]:
        """
        Extract the Murcko scaffold (core) from a molecule.

        Args:
            mol_smiles: SMILES of the molecule

        Returns:
            SMILES of the core scaffold, or None if extraction fails
        """
        mol = Chem.MolFromSmiles(mol_smiles)
        if mol is None:
            return None

        try:
            core = MurckoScaffold.GetScaffoldForMol(mol)
            return Chem.MolToSmiles(core)
        except Exception:
            return None

    def extract_generic_core(self, mol_smiles: str) -> Optional[str]:
        """
        Extract the generic Murcko scaffold (all atoms replaced with carbons).

        Args:
            mol_smiles: SMILES of the molecule

        Returns:
            SMILES of the generic core scaffold
        """
        mol = Chem.MolFromSmiles(mol_smiles)
        if mol is None:
            return None

        try:
            generic_core = MurckoScaffold.MakeScaffoldGeneric(
                MurckoScaffold.GetScaffoldForMol(mol)
            )
            return Chem.MolToSmiles(generic_core)
        except Exception:
            return None

    def decompose_molecule(
        self, mol_smiles: str
    ) -> Optional[Dict[str, Any]]:
        """
        Decompose a molecule into core and side chains (R-groups).

        Args:
            mol_smiles: SMILES of the molecule to decompose

        Returns:
            Dictionary with 'core' and 'rgroups' keys
        """
        mol = Chem.MolFromSmiles(mol_smiles)
        if mol is None:
            return None

        try:
            core_mol = MurckoScaffold.GetScaffoldForMol(mol)
            core_smiles = Chem.MolToSmiles(core_mol)

            # Find atoms in the core
            core_match = mol.GetSubstructMatch(core_mol)

            # Extract R-groups (atoms not in core)
            rgroup_atoms = set(range(mol.GetNumAtoms())) - set(core_match)

            return {
                "core": core_smiles,
                "core_mol": core_mol,
                "core_atoms": list(core_match),
                "rgroup_atoms": list(rgroup_atoms),
            }
        except Exception:
            return None

    def find_hopping_candidates(
        self,
        mol_smiles: str,
        similarity_threshold: float = 0.3,
        max_candidates: int = 20,
        require_ring_count_match: bool = False
    ) -> List[CoreCandidate]:
        """
        Find candidate cores for scaffold hopping.

        Args:
            mol_smiles: SMILES of the query molecule
            similarity_threshold: Minimum similarity for candidates
            max_candidates: Maximum number of candidates to return
            require_ring_count_match: If True, only return cores with same ring count

        Returns:
            List of CoreCandidate objects
        """
        # Extract core from query molecule
        core_smiles = self.extract_core(mol_smiles)
        if core_smiles is None:
            return []

        # Find similar cores in database
        similar_cores = self.database.find_similar(
            core_smiles,
            threshold=similarity_threshold,
            top_n=max_candidates * 2
        )

        candidates = []
        query_mol = Chem.MolFromSmiles(core_smiles)
        query_ring_count = query_mol.GetRingInfo().NumRings() if query_mol else 0

        for smiles, similarity in similar_cores:
            mol = self.database.lookup(smiles)
            if mol is None:
                continue

            # Check ring count if required
            if require_ring_count_match:
                ring_count = mol.GetRingInfo().NumRings()
                if ring_count != query_ring_count:
                    continue

            # Calculate properties
            properties = self._calculate_core_properties(mol)

            # Find potential attachment points
            attachment_points = self._find_attachment_points(mol)

            candidate = CoreCandidate(
                smiles=smiles,
                mol=mol,
                similarity_score=similarity,
                attachment_points=attachment_points,
                properties=properties,
                source="database"
            )
            candidates.append(candidate)

            if len(candidates) >= max_candidates:
                break

        return candidates

    def _calculate_core_properties(self, mol: Any) -> Dict[str, float]:
        """Calculate properties for a core scaffold."""
        return {
            "mw": Descriptors.MolWt(mol),
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "rings": mol.GetRingInfo().NumRings(),
            "aromatic_rings": Descriptors.NumAromaticRings(mol),
            "heteroatoms": Descriptors.NumHeteroatoms(mol),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        }

    def _find_attachment_points(self, mol: Any) -> List[int]:
        """
        Find potential attachment points on a core scaffold.

        These are atoms that could serve as connection points for R-groups.
        Typically includes atoms with available valence for substitution.
        """
        attachment_points = []

        for atom in mol.GetAtoms():
            # Check for dummy atoms (explicit attachment points)
            if atom.GetAtomicNum() == 0:
                attachment_points.append(atom.GetIdx())
                continue

            # Check for atoms with available valence
            implicit_h = atom.GetTotalNumHs()
            if implicit_h > 0:
                # This atom has hydrogens that could be replaced
                attachment_points.append(atom.GetIdx())

        return attachment_points

    def get_compatible_cores(
        self,
        mol_smiles: str,
        required_attachment_count: int
    ) -> List[CoreCandidate]:
        """
        Find cores that have the required number of attachment points.

        Args:
            mol_smiles: SMILES of the query molecule
            required_attachment_count: Number of attachment points needed

        Returns:
            List of CoreCandidate objects with matching attachment count
        """
        candidates = self.find_hopping_candidates(mol_smiles, max_candidates=50)

        compatible = [
            c for c in candidates
            if len(c.attachment_points) >= required_attachment_count
        ]

        return compatible

    def calculate_hopping_potential(
        self,
        original_core: str,
        candidate_core: str
    ) -> Dict[str, float]:
        """
        Calculate metrics indicating the hopping potential between two cores.

        Args:
            original_core: SMILES of the original core
            candidate_core: SMILES of the candidate core

        Returns:
            Dictionary of hopping potential metrics
        """
        orig_mol = Chem.MolFromSmiles(original_core)
        cand_mol = Chem.MolFromSmiles(candidate_core)

        if orig_mol is None or cand_mol is None:
            return {}

        # Calculate fingerprint similarity
        orig_fp = FingerprintMols.FingerprintMol(orig_mol)
        cand_fp = FingerprintMols.FingerprintMol(cand_mol)
        similarity = DataStructs.TanimotoSimilarity(orig_fp, cand_fp)

        # Calculate property differences
        orig_props = self._calculate_core_properties(orig_mol)
        cand_props = self._calculate_core_properties(cand_mol)

        return {
            "similarity": similarity,
            "novelty": 1.0 - similarity,  # Higher novelty = more different
            "mw_diff": abs(orig_props["mw"] - cand_props["mw"]),
            "ring_diff": abs(orig_props["rings"] - cand_props["rings"]),
            "size_diff": abs(orig_props["heavy_atoms"] - cand_props["heavy_atoms"]),
        }
