"""Molecular similarity calculation utilities."""

from typing import Optional

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdMolDescriptors


class SimilarityCalculator:
    """Calculate various types of molecular similarity."""

    def __init__(self, fingerprint_type: str = "morgan", radius: int = 2, n_bits: int = 2048):
        """
        Initialize the similarity calculator.

        Args:
            fingerprint_type: Type of fingerprint ('morgan', 'rdkit', 'maccs')
            radius: Radius for Morgan fingerprints
            n_bits: Number of bits for fingerprint
        """
        self.fingerprint_type = fingerprint_type
        self.radius = radius
        self.n_bits = n_bits

    def get_fingerprint(self, mol: Chem.Mol) -> Optional[DataStructs.ExplicitBitVect]:
        """
        Generate fingerprint for a molecule.

        Args:
            mol: RDKit molecule object

        Returns:
            Fingerprint bit vector or None if failed
        """
        if mol is None:
            return None

        if self.fingerprint_type == "morgan":
            return AllChem.GetMorganFingerprintAsBitVect(mol, self.radius, nBits=self.n_bits)
        elif self.fingerprint_type == "rdkit":
            return Chem.RDKFingerprint(mol, fpSize=self.n_bits)
        elif self.fingerprint_type == "maccs":
            return rdMolDescriptors.GetMACCSKeysFingerprint(mol)
        else:
            raise ValueError(f"Unknown fingerprint type: {self.fingerprint_type}")

    def tanimoto_similarity(self, mol1: Chem.Mol, mol2: Chem.Mol) -> float:
        """
        Calculate Tanimoto similarity between two molecules.

        Args:
            mol1: First molecule
            mol2: Second molecule

        Returns:
            Tanimoto similarity score (0-1)
        """
        fp1 = self.get_fingerprint(mol1)
        fp2 = self.get_fingerprint(mol2)

        if fp1 is None or fp2 is None:
            return 0.0

        return DataStructs.TanimotoSimilarity(fp1, fp2)

    def shape_similarity(self, mol1: Chem.Mol, mol2: Chem.Mol) -> float:
        """
        Calculate shape similarity between two molecules using 3D conformers.

        Args:
            mol1: First molecule (with or without 3D coords)
            mol2: Second molecule (with or without 3D coords)

        Returns:
            Shape similarity score (0-1)
        """
        # Generate 3D conformers if needed
        mol1_3d = self._ensure_3d(mol1)
        mol2_3d = self._ensure_3d(mol2)

        if mol1_3d is None or mol2_3d is None:
            return 0.0

        try:
            from rdkit.Chem import AllChem, rdShapeHelpers

            # Align molecules and calculate shape similarity
            align_score = AllChem.GetO3A(mol1_3d, mol2_3d).Align()
            shape_sim = rdShapeHelpers.ShapeTanimotoDist(mol1_3d, mol2_3d)
            # Convert distance to similarity
            return 1.0 - shape_sim
        except Exception:
            # Fallback if shape calculation fails
            return 0.0

    def _ensure_3d(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """
        Ensure molecule has 3D coordinates.

        Args:
            mol: RDKit molecule

        Returns:
            Molecule with 3D coordinates or None if failed
        """
        if mol is None:
            return None

        mol_copy = Chem.Mol(mol)

        try:
            # Check if already has 3D coords
            conf = mol_copy.GetConformer() if mol_copy.GetNumConformers() > 0 else None
            if conf is None or not conf.Is3D():
                mol_copy = Chem.AddHs(mol_copy)
                AllChem.EmbedMolecule(mol_copy, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(mol_copy)
            return mol_copy
        except Exception:
            return None

    def esp_similarity(self, mol1: Chem.Mol, mol2: Chem.Mol) -> float:
        """
        Calculate Electrostatic Potential (ESP) similarity between molecules.

        Uses partial charge distribution as a proxy for ESP.

        Args:
            mol1: First molecule
            mol2: Second molecule

        Returns:
            ESP similarity score (0-1)
        """
        if mol1 is None or mol2 is None:
            return 0.0

        try:
            # Compute Gasteiger charges as ESP proxy
            mol1_copy = Chem.Mol(mol1)
            mol2_copy = Chem.Mol(mol2)

            AllChem.ComputeGasteigerCharges(mol1_copy)
            AllChem.ComputeGasteigerCharges(mol2_copy)

            # Get charge distributions
            charges1 = []
            for atom in mol1_copy.GetAtoms():
                charge = atom.GetDoubleProp("_GasteigerCharge")
                if not (charge != charge):  # Check for NaN
                    charges1.append(charge)

            charges2 = []
            for atom in mol2_copy.GetAtoms():
                charge = atom.GetDoubleProp("_GasteigerCharge")
                if not (charge != charge):  # Check for NaN
                    charges2.append(charge)

            if not charges1 or not charges2:
                return 0.0

            # Calculate statistics for comparison
            avg1 = sum(charges1) / len(charges1)
            avg2 = sum(charges2) / len(charges2)

            std1 = (sum((c - avg1) ** 2 for c in charges1) / len(charges1)) ** 0.5
            std2 = (sum((c - avg2) ** 2 for c in charges2) / len(charges2)) ** 0.5

            # Compare charge distributions
            avg_diff = abs(avg1 - avg2)
            std_diff = abs(std1 - std2) if std1 > 0 and std2 > 0 else 0

            # Normalize to 0-1 similarity
            similarity = max(0.0, 1.0 - (avg_diff + std_diff) / 2)
            return similarity

        except Exception:
            return 0.0

    def combined_similarity(
        self,
        mol1: Chem.Mol,
        mol2: Chem.Mol,
        fingerprint_weight: float = 0.5,
        shape_weight: float = 0.3,
        esp_weight: float = 0.2,
    ) -> float:
        """
        Calculate combined similarity using multiple metrics.

        Args:
            mol1: First molecule
            mol2: Second molecule
            fingerprint_weight: Weight for fingerprint similarity
            shape_weight: Weight for shape similarity
            esp_weight: Weight for ESP similarity

        Returns:
            Combined similarity score (0-1)
        """
        fp_sim = self.tanimoto_similarity(mol1, mol2)
        shape_sim = self.shape_similarity(mol1, mol2)
        esp_sim = self.esp_similarity(mol1, mol2)

        total_weight = fingerprint_weight + shape_weight + esp_weight
        return (
            fingerprint_weight * fp_sim + shape_weight * shape_sim + esp_weight * esp_sim
        ) / total_weight

    def rank_by_similarity(
        self,
        candidates: list[Chem.Mol],
        target_mol: Chem.Mol,
        method: str = "tanimoto",
        top_n: Optional[int] = None,
    ) -> list[tuple[Chem.Mol, float]]:
        """
        Rank candidates by similarity to target molecule.

        Args:
            candidates: List of candidate molecules
            target_mol: Target molecule
            method: Similarity method ('tanimoto', 'shape', 'esp', 'combined')
            top_n: Return only top N results

        Returns:
            List of (molecule, score) tuples sorted by score descending
        """
        scored = []

        for mol in candidates:
            if mol is None:
                continue

            if method == "tanimoto":
                score = self.tanimoto_similarity(mol, target_mol)
            elif method == "shape":
                score = self.shape_similarity(mol, target_mol)
            elif method == "esp":
                score = self.esp_similarity(mol, target_mol)
            elif method == "combined":
                score = self.combined_similarity(mol, target_mol)
            else:
                raise ValueError(f"Unknown similarity method: {method}")

            scored.append((mol, score))

        scored.sort(key=lambda x: x[1], reverse=True)

        if top_n is not None:
            scored = scored[:top_n]

        return scored
