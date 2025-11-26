"""CoreHopping workflow implementation."""

from dataclasses import dataclass
from typing import Optional

from rdkit import Chem

from expopt.database.interface import DatabaseInterface
from expopt.database.mock_database import MockDatabase
from expopt.utils.constraints import AtomCountConstraint, ConnectionSiteConstraint
from expopt.utils.scoring import PhysicochemicalScorer
from expopt.utils.similarity import SimilarityCalculator


@dataclass
class CoreCandidate:
    """Represents a core hopping candidate with scores."""

    smiles: str
    mol: Chem.Mol
    connection_sites: int
    physicochemical_score: float
    shape_score: Optional[float] = None
    esp_score: Optional[float] = None
    combined_score: Optional[float] = None
    source: str = "database"


class CoreHopping:
    """
    CoreHopping workflow for replacing core skeletons.

    Maintains side chains and key functional groups while
    finding alternative cores with hopping potential.
    """

    def __init__(
        self,
        database: Optional[DatabaseInterface] = None,
        top_n_physicochemical: int = 1000,
        top_n_final: int = 100,
    ):
        """
        Initialize CoreHopping workflow.

        Args:
            database: Database interface for core lookup
            top_n_physicochemical: Number of candidates after physicochemical filtering
            top_n_final: Number of final candidates
        """
        self.database = database or MockDatabase()
        self.top_n_physicochemical = top_n_physicochemical
        self.top_n_final = top_n_final
        self.scorer = PhysicochemicalScorer()
        self.similarity = SimilarityCalculator()

    def find_alternative_cores(
        self,
        target_core: str,
        min_connection_sites: Optional[int] = None,
        use_fallback: bool = True,
    ) -> list[CoreCandidate]:
        """
        Find alternative cores for hopping.

        Args:
            target_core: SMILES of target core to replace
            min_connection_sites: Minimum required connection sites (default: target sites + 1)
            use_fallback: Whether to use fallback generation

        Returns:
            List of CoreCandidate objects sorted by score
        """
        target_mol = Chem.MolFromSmiles(target_core)
        if target_mol is None:
            raise ValueError(f"Invalid target core SMILES: {target_core}")

        # Determine minimum connection sites
        site_constraint = ConnectionSiteConstraint(1)
        target_sites = site_constraint.count_connection_sites(target_mol)
        min_sites = min_connection_sites or target_sites

        # Stage 1: Database search
        candidates = self._database_search(target_core, target_mol, min_sites)

        # Stage 2: Fallback if needed
        if not candidates and use_fallback:
            candidates = self._fallback_generation(target_mol, min_sites)

        # Apply constraints and rank
        candidates = self._apply_constraints(candidates, target_mol, min_sites)
        return self._rank_candidates(candidates, target_mol)

    def _database_search(
        self,
        target_core: str,
        target_mol: Chem.Mol,
        min_sites: int,
    ) -> list[CoreCandidate]:
        """
        Search database for alternative cores.

        Args:
            target_core: SMILES of target core
            target_mol: Target core as RDKit Mol
            min_sites: Minimum connection sites required

        Returns:
            List of candidates from database
        """
        candidates = []
        seen_smiles = {target_core}

        # Get alternative cores from database
        alt_cores = self.database.get_alternative_cores(target_core)

        for core_record in alt_cores:
            if core_record.core_smiles not in seen_smiles:
                mol = Chem.MolFromSmiles(core_record.core_smiles)
                if mol is not None and core_record.connection_sites >= min_sites:
                    score = self.scorer.calculate_score(mol, target_mol)
                    candidates.append(
                        CoreCandidate(
                            smiles=core_record.core_smiles,
                            mol=mol,
                            connection_sites=core_record.connection_sites,
                            physicochemical_score=score,
                            source="database",
                        )
                    )
                    seen_smiles.add(core_record.core_smiles)

        # Also search cores with minimum sites
        cores = self.database.search_cores(min_sites)
        for core_record in cores:
            if core_record.core_smiles not in seen_smiles:
                mol = Chem.MolFromSmiles(core_record.core_smiles)
                if mol is not None:
                    score = self.scorer.calculate_score(mol, target_mol)
                    candidates.append(
                        CoreCandidate(
                            smiles=core_record.core_smiles,
                            mol=mol,
                            connection_sites=core_record.connection_sites,
                            physicochemical_score=score,
                            source="database",
                        )
                    )
                    seen_smiles.add(core_record.core_smiles)

        return candidates

    def _fallback_generation(
        self,
        target_mol: Chem.Mol,
        min_sites: int,
    ) -> list[CoreCandidate]:
        """
        Generate candidates when no database hits.

        Args:
            target_mol: Target core
            min_sites: Minimum connection sites

        Returns:
            List of generated candidates
        """
        candidates = []

        # Common core scaffolds for fallback
        common_cores = [
            ("c1ccccc1", 6),  # Benzene
            ("c1ccncc1", 5),  # Pyridine
            ("c1ccc2c(c1)cccc2", 8),  # Naphthalene
            ("c1ccc2[nH]ccc2c1", 7),  # Indole
            ("c1ccc2ncccc2c1", 7),  # Quinoline
            ("c1ccc2occc2c1", 6),  # Benzofuran
            ("c1ccc2sccc2c1", 6),  # Benzothiophene
            ("c1cnc2ccccc2n1", 6),  # Quinazoline
            ("c1ccc2c(c1)nccn2", 6),  # Quinoxaline
            ("c1ccc2c(c1)c3ccccc3cc2", 10),  # Anthracene
            ("c1ccoc1", 4),  # Furan
            ("c1ccsc1", 4),  # Thiophene
            ("c1cc[nH]c1", 4),  # Pyrrole
            ("c1cnc[nH]1", 3),  # Imidazole
            ("c1ccnc(n1)N", 4),  # 2-aminopyrimidine
        ]

        # Apply atom count constraint
        target_atom_count = target_mol.GetNumHeavyAtoms()
        atom_constraint = AtomCountConstraint(target_atom_count)
        site_constraint = ConnectionSiteConstraint(min_sites)

        for core_smiles, sites in common_cores:
            mol = Chem.MolFromSmiles(core_smiles)
            if mol is not None:
                if atom_constraint.check_constraint(mol) and sites >= min_sites:
                    score = self.scorer.calculate_score(mol, target_mol)
                    candidates.append(
                        CoreCandidate(
                            smiles=core_smiles,
                            mol=mol,
                            connection_sites=sites,
                            physicochemical_score=score,
                            source="generated",
                        )
                    )

        # Sort by score and keep top N
        candidates.sort(key=lambda c: c.physicochemical_score, reverse=True)
        return candidates[: self.top_n_physicochemical]

    def _apply_constraints(
        self,
        candidates: list[CoreCandidate],
        target_mol: Chem.Mol,
        min_sites: int,
    ) -> list[CoreCandidate]:
        """
        Apply atom count and connection site constraints.

        Args:
            candidates: List of candidates
            target_mol: Target core
            min_sites: Minimum connection sites

        Returns:
            Filtered list of candidates
        """
        target_atom_count = target_mol.GetNumHeavyAtoms()
        atom_constraint = AtomCountConstraint(target_atom_count)

        filtered = []
        for candidate in candidates:
            if (
                atom_constraint.check_constraint(candidate.mol)
                and candidate.connection_sites >= min_sites
            ):
                filtered.append(candidate)

        return filtered

    def _rank_candidates(
        self,
        candidates: list[CoreCandidate],
        target_mol: Chem.Mol,
    ) -> list[CoreCandidate]:
        """
        Rank candidates using shape and ESP similarity.

        Args:
            candidates: List of candidates
            target_mol: Target core

        Returns:
            Ranked list of candidates
        """
        if not candidates:
            return []

        # Calculate shape and ESP scores
        for candidate in candidates[: self.top_n_physicochemical]:
            candidate.shape_score = self.similarity.shape_similarity(
                candidate.mol, target_mol
            )
            candidate.esp_score = self.similarity.esp_similarity(
                candidate.mol, target_mol
            )

            # Combined score
            candidate.combined_score = (
                0.4 * candidate.physicochemical_score
                + 0.35 * (candidate.shape_score or 0)
                + 0.25 * (candidate.esp_score or 0)
            )

        # Sort by combined score
        candidates.sort(key=lambda c: c.combined_score or 0, reverse=True)
        return candidates[: self.top_n_final]

    def get_alternative_core_molecules(
        self,
        target_core: str,
        min_connection_sites: Optional[int] = None,
    ) -> list[tuple[str, float]]:
        """
        Convenience method to get alternative core SMILES with scores.

        Args:
            target_core: SMILES of target core
            min_connection_sites: Minimum connection sites required

        Returns:
            List of (SMILES, score) tuples
        """
        candidates = self.find_alternative_cores(target_core, min_connection_sites)
        return [
            (c.smiles, c.combined_score or c.physicochemical_score) for c in candidates
        ]
