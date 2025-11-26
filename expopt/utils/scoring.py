"""Physicochemical property scoring utilities."""

from typing import Optional

from rdkit import Chem
from rdkit.Chem import Descriptors


class PhysicochemicalScorer:
    """Calculate and score physicochemical properties of molecules/fragments."""

    def __init__(self):
        """Initialize the scorer."""
        pass

    def calculate_properties(self, mol: Chem.Mol) -> dict:
        """
        Calculate physicochemical properties for a molecule.

        Args:
            mol: RDKit molecule object

        Returns:
            Dictionary containing calculated properties
        """
        if mol is None:
            return {}

        return {
            "mw": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "tpsa": Descriptors.TPSA(mol),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "num_atoms": mol.GetNumAtoms(),
            "num_heavy_atoms": mol.GetNumHeavyAtoms(),
            "num_rings": Descriptors.RingCount(mol),
            "num_aromatic_rings": Descriptors.NumAromaticRings(mol),
        }

    def calculate_score(
        self,
        candidate_mol: Chem.Mol,
        target_mol: Chem.Mol,
        weights: Optional[dict] = None,
    ) -> float:
        """
        Calculate similarity score between candidate and target based on properties.

        Args:
            candidate_mol: Candidate molecule
            target_mol: Target molecule to compare against
            weights: Optional dictionary of property weights

        Returns:
            Float score representing property similarity (0-1, higher is better)
        """
        if candidate_mol is None or target_mol is None:
            return 0.0

        default_weights = {
            "mw": 0.2,
            "logp": 0.25,
            "hbd": 0.1,
            "hba": 0.1,
            "tpsa": 0.15,
            "rotatable_bonds": 0.05,
            "num_heavy_atoms": 0.1,
            "num_aromatic_rings": 0.05,
        }

        weights = weights or default_weights

        candidate_props = self.calculate_properties(candidate_mol)
        target_props = self.calculate_properties(target_mol)

        score = 0.0
        total_weight = 0.0

        for prop, weight in weights.items():
            if prop in candidate_props and prop in target_props:
                target_val = target_props[prop]
                candidate_val = candidate_props[prop]

                if target_val == 0:
                    prop_score = 1.0 if candidate_val == 0 else 0.0
                else:
                    # Calculate similarity based on relative difference
                    diff = abs(candidate_val - target_val)
                    prop_score = max(0.0, 1.0 - diff / (abs(target_val) + 1))

                score += weight * prop_score
                total_weight += weight

        return score / total_weight if total_weight > 0 else 0.0

    def rank_candidates(
        self,
        candidates: list[Chem.Mol],
        target_mol: Chem.Mol,
        weights: Optional[dict] = None,
        top_n: Optional[int] = None,
    ) -> list[tuple[Chem.Mol, float]]:
        """
        Rank candidate molecules by physicochemical property similarity.

        Args:
            candidates: List of candidate molecules
            target_mol: Target molecule to compare against
            weights: Optional property weights
            top_n: Return only top N results

        Returns:
            List of (molecule, score) tuples sorted by score descending
        """
        scored = []
        for mol in candidates:
            if mol is not None:
                score = self.calculate_score(mol, target_mol, weights)
                scored.append((mol, score))

        scored.sort(key=lambda x: x[1], reverse=True)

        if top_n is not None:
            scored = scored[:top_n]

        return scored
