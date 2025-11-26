"""
RGroup Replacement Module

This module provides functionality for finding candidate fragments to replace
target fragments while ensuring similar properties and ligand site matching.

Workflow:
1. Priority: Exact matching (database hit)
2. Fallback: Database miss (candidate generation and filtering)
"""

from typing import List, Optional, Dict, Any, Tuple
from dataclasses import dataclass

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdFMCS
    from rdkit.Chem.Fingerprints import FingerprintMols
    from rdkit import DataStructs
except ImportError:
    Chem = None
    AllChem = None
    Descriptors = None
    rdFMCS = None
    FingerprintMols = None
    DataStructs = None


@dataclass
class FragmentCandidate:
    """Represents a candidate fragment for replacement."""
    smiles: str
    mol: Any  # RDKit Mol object
    similarity_score: float = 0.0
    properties: Optional[Dict[str, float]] = None
    source: str = "generated"  # "database" or "generated"

    def __post_init__(self):
        if self.properties is None:
            self.properties = {}


class FragmentDatabase:
    """Database of fragments for exact matching."""

    def __init__(self):
        self._fragments: Dict[str, Any] = {}  # canonical SMILES -> Mol object
        self._fingerprints: Dict[str, Any] = {}

    def add_fragment(self, smiles: str) -> bool:
        """Add a fragment to the database."""
        if Chem is None:
            raise ImportError("RDKit is required for this functionality")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        canonical_smiles = Chem.MolToSmiles(mol)
        self._fragments[canonical_smiles] = mol
        self._fingerprints[canonical_smiles] = FingerprintMols.FingerprintMol(mol)
        return True

    def add_fragments_from_list(self, smiles_list: List[str]) -> int:
        """Add multiple fragments to the database. Returns count of successfully added."""
        count = 0
        for smiles in smiles_list:
            if self.add_fragment(smiles):
                count += 1
        return count

    def lookup(self, smiles: str) -> Optional[Any]:
        """Look up a fragment by SMILES. Returns Mol object if found."""
        if Chem is None:
            raise ImportError("RDKit is required for this functionality")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        canonical_smiles = Chem.MolToSmiles(mol)
        return self._fragments.get(canonical_smiles)

    def find_similar(
        self, smiles: str, threshold: float = 0.5, top_n: int = 10
    ) -> List[Tuple[str, float]]:
        """Find fragments similar to the query. Returns list of (smiles, similarity)."""
        if Chem is None:
            raise ImportError("RDKit is required for this functionality")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []

        query_fp = FingerprintMols.FingerprintMol(mol)
        results = []

        for frag_smiles, frag_fp in self._fingerprints.items():
            similarity = DataStructs.TanimotoSimilarity(query_fp, frag_fp)
            if similarity >= threshold:
                results.append((frag_smiles, similarity))

        results.sort(key=lambda x: x[1], reverse=True)
        return results[:top_n]

    def get_all_fragments(self) -> List[str]:
        """Return all fragment SMILES in the database."""
        return list(self._fragments.keys())

    def __len__(self) -> int:
        return len(self._fragments)


class RGroupReplacer:
    """
    RGroup Replacement Engine

    Finds candidate fragments to replace a target fragment while maintaining
    similar properties and ligand site matching.
    """

    def __init__(self, database: Optional[FragmentDatabase] = None):
        """
        Initialize the RGroup replacer.

        Args:
            database: Optional fragment database for exact matching
        """
        if Chem is None:
            raise ImportError("RDKit is required for this functionality")

        self.database = database or FragmentDatabase()

    def exact_match(self, target_smiles: str) -> Optional[FragmentCandidate]:
        """
        Attempt exact match in the database.

        Args:
            target_smiles: SMILES of the target fragment to match

        Returns:
            FragmentCandidate if exact match found, None otherwise
        """
        mol = self.database.lookup(target_smiles)
        if mol is not None:
            return FragmentCandidate(
                smiles=Chem.MolToSmiles(mol),
                mol=mol,
                similarity_score=1.0,
                source="database"
            )
        return None

    def generate_candidates(
        self,
        target_smiles: str,
        similarity_threshold: float = 0.5,
        max_candidates: int = 10,
        property_constraints: Optional[Dict[str, Tuple[float, float]]] = None
    ) -> List[FragmentCandidate]:
        """
        Generate candidate fragments when exact match is not found.

        This method finds similar fragments from the database and filters
        them based on property constraints.

        Args:
            target_smiles: SMILES of the target fragment
            similarity_threshold: Minimum Tanimoto similarity score
            max_candidates: Maximum number of candidates to return
            property_constraints: Dict of property name -> (min, max) constraints

        Returns:
            List of FragmentCandidate objects
        """
        similar_fragments = self.database.find_similar(
            target_smiles,
            threshold=similarity_threshold,
            top_n=max_candidates * 2  # Get more to filter
        )

        candidates = []
        for smiles, similarity in similar_fragments:
            mol = self.database.lookup(smiles)
            if mol is None:
                continue

            # Calculate properties
            properties = self._calculate_properties(mol)

            # Apply property constraints if specified
            if property_constraints is not None:
                if not self._check_property_constraints(properties, property_constraints):
                    continue

            candidate = FragmentCandidate(
                smiles=smiles,
                mol=mol,
                similarity_score=similarity,
                properties=properties,
                source="generated"
            )
            candidates.append(candidate)

            if len(candidates) >= max_candidates:
                break

        return candidates

    def find_replacements(
        self,
        target_smiles: str,
        similarity_threshold: float = 0.5,
        max_candidates: int = 10,
        property_constraints: Optional[Dict[str, Tuple[float, float]]] = None
    ) -> List[FragmentCandidate]:
        """
        Find replacement fragments for the target.

        Workflow:
        1. Try exact match first
        2. Fall back to candidate generation if no exact match

        Args:
            target_smiles: SMILES of the target fragment
            similarity_threshold: Minimum similarity for fallback search
            max_candidates: Maximum number of candidates to return
            property_constraints: Property constraints for filtering

        Returns:
            List of FragmentCandidate objects, exact match first if found
        """
        candidates = []

        # Priority: Try exact match
        exact = self.exact_match(target_smiles)
        if exact is not None:
            candidates.append(exact)

        # Fallback: Generate similar candidates
        generated = self.generate_candidates(
            target_smiles,
            similarity_threshold=similarity_threshold,
            max_candidates=max_candidates - len(candidates),
            property_constraints=property_constraints
        )
        candidates.extend(generated)

        return candidates

    def _calculate_properties(self, mol: Any) -> Dict[str, float]:
        """Calculate molecular properties for a fragment."""
        return {
            "mw": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "tpsa": Descriptors.TPSA(mol),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "heavy_atoms": mol.GetNumHeavyAtoms(),
        }

    def _check_property_constraints(
        self,
        properties: Dict[str, float],
        constraints: Dict[str, Tuple[float, float]]
    ) -> bool:
        """Check if properties satisfy constraints."""
        for prop_name, (min_val, max_val) in constraints.items():
            if prop_name in properties:
                value = properties[prop_name]
                if value < min_val or value > max_val:
                    return False
        return True

    def get_attachment_points(self, smiles: str) -> List[int]:
        """
        Get the indices of atoms that serve as attachment points.

        Attachment points are typically indicated by dummy atoms (*)
        or atoms with explicit attachment markers.

        Args:
            smiles: SMILES string of the fragment

        Returns:
            List of atom indices that are attachment points
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []

        attachment_points = []
        for atom in mol.GetAtoms():
            # Check for dummy atoms (atomic number 0) or isotope labels
            if atom.GetAtomicNum() == 0:
                # Find the real atom connected to the dummy
                for neighbor in atom.GetNeighbors():
                    attachment_points.append(neighbor.GetIdx())
            # Also check for common attachment point markers (R groups)
            elif atom.GetSymbol() == "*":
                for neighbor in atom.GetNeighbors():
                    attachment_points.append(neighbor.GetIdx())

        return list(set(attachment_points))
