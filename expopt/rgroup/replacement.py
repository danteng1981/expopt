"""RGroup replacement workflow implementation."""

from dataclasses import dataclass
from typing import Optional

from rdkit import Chem

from expopt.database.interface import DatabaseInterface
from expopt.database.mock_database import MockDatabase
from expopt.utils.constraints import AtomCountConstraint
from expopt.utils.scoring import PhysicochemicalScorer
from expopt.utils.similarity import SimilarityCalculator


@dataclass
class ReplacementCandidate:
    """Represents a replacement candidate with scores."""

    smiles: str
    mol: Chem.Mol
    physicochemical_score: float
    shape_score: Optional[float] = None
    esp_score: Optional[float] = None
    combined_score: Optional[float] = None
    source: str = "database"  # 'database' or 'generated'


class RGroupReplacement:
    """
    RGroup replacement workflow for finding candidate fragments.

    Implements two-stage workflow:
    1. Precise Matching: Search database for patent records and scaffolds
    2. Fallback: Generate and filter candidates based on constraints and scoring
    """

    def __init__(
        self,
        database: Optional[DatabaseInterface] = None,
        top_n_physicochemical: int = 1000,
        top_n_final: int = 100,
    ):
        """
        Initialize RGroup replacement workflow.

        Args:
            database: Database interface for patent/scaffold lookup
            top_n_physicochemical: Number of candidates to keep after physicochemical scoring
            top_n_final: Number of final candidates after shape/ESP refinement
        """
        self.database = database or MockDatabase()
        self.top_n_physicochemical = top_n_physicochemical
        self.top_n_final = top_n_final
        self.scorer = PhysicochemicalScorer()
        self.similarity = SimilarityCalculator()

    def find_replacements(
        self,
        target_fragment: str,
        scaffold_smiles: Optional[str] = None,
        use_fallback: bool = True,
    ) -> list[ReplacementCandidate]:
        """
        Find replacement candidates for a target fragment.

        Args:
            target_fragment: SMILES of target fragment to replace
            scaffold_smiles: Optional SMILES of scaffold context
            use_fallback: Whether to use fallback generation if no database hits

        Returns:
            List of ReplacementCandidate objects sorted by score
        """
        target_mol = Chem.MolFromSmiles(target_fragment)
        if target_mol is None:
            raise ValueError(f"Invalid target fragment SMILES: {target_fragment}")

        # Stage 1: Precise Matching (Database search)
        candidates = self._database_search(target_fragment, scaffold_smiles, target_mol)

        # Stage 2: Fallback if no database hits
        if not candidates and use_fallback:
            candidates = self._fallback_generation(target_mol)

        # Final ranking
        return self._rank_candidates(candidates, target_mol)

    def _database_search(
        self,
        target_fragment: str,
        scaffold_smiles: Optional[str],
        target_mol: Chem.Mol,
    ) -> list[ReplacementCandidate]:
        """
        Search database for replacement candidates.

        Args:
            target_fragment: SMILES of target fragment
            scaffold_smiles: Optional scaffold context
            target_mol: Target fragment as RDKit Mol

        Returns:
            List of candidates from database
        """
        candidates = []
        seen_smiles = {target_fragment}  # Don't include target itself

        # Search by scaffold if provided
        if scaffold_smiles:
            alternative_fragments = self.database.get_alternative_fragments(scaffold_smiles)
            for frag_smiles in alternative_fragments:
                if frag_smiles not in seen_smiles:
                    mol = Chem.MolFromSmiles(frag_smiles)
                    if mol is not None:
                        score = self.scorer.calculate_score(mol, target_mol)
                        candidates.append(
                            ReplacementCandidate(
                                smiles=frag_smiles,
                                mol=mol,
                                physicochemical_score=score,
                                source="database",
                            )
                        )
                        seen_smiles.add(frag_smiles)

        # Search by fragment in patents
        patents = self.database.search_by_fragment(target_fragment)
        for patent in patents:
            if patent.fragments:
                for frag_smiles in patent.fragments:
                    if frag_smiles not in seen_smiles:
                        mol = Chem.MolFromSmiles(frag_smiles)
                        if mol is not None:
                            score = self.scorer.calculate_score(mol, target_mol)
                            candidates.append(
                                ReplacementCandidate(
                                    smiles=frag_smiles,
                                    mol=mol,
                                    physicochemical_score=score,
                                    source="database",
                                )
                            )
                            seen_smiles.add(frag_smiles)

        return candidates

    def _fallback_generation(self, target_mol: Chem.Mol) -> list[ReplacementCandidate]:
        """
        Generate candidates when no database hits found.

        Args:
            target_mol: Target fragment as RDKit Mol

        Returns:
            List of generated candidates
        """
        candidates = []

        # Get fragment library from database
        if hasattr(self.database, "get_fragment_library"):
            fragment_library = self.database.get_fragment_library()
        else:
            fragment_library = []

        # Apply atom count constraint
        target_atom_count = target_mol.GetNumHeavyAtoms()
        constraint = AtomCountConstraint(target_atom_count)

        for frag_smiles in fragment_library:
            mol = Chem.MolFromSmiles(frag_smiles)
            if mol is not None and constraint.check_constraint(mol):
                score = self.scorer.calculate_score(mol, target_mol)
                candidates.append(
                    ReplacementCandidate(
                        smiles=frag_smiles,
                        mol=mol,
                        physicochemical_score=score,
                        source="generated",
                    )
                )

        # Sort by physicochemical score and keep top N
        candidates.sort(key=lambda c: c.physicochemical_score, reverse=True)
        return candidates[: self.top_n_physicochemical]

    def _rank_candidates(
        self,
        candidates: list[ReplacementCandidate],
        target_mol: Chem.Mol,
    ) -> list[ReplacementCandidate]:
        """
        Rank candidates using shape and ESP similarity.

        Args:
            candidates: List of candidates to rank
            target_mol: Target fragment

        Returns:
            Ranked list of candidates
        """
        if not candidates:
            return []

        # Calculate shape and ESP scores for top candidates
        for candidate in candidates[: self.top_n_physicochemical]:
            candidate.shape_score = self.similarity.shape_similarity(
                candidate.mol, target_mol
            )
            candidate.esp_score = self.similarity.esp_similarity(candidate.mol, target_mol)

            # Combined score (weighted average)
            candidate.combined_score = (
                0.4 * candidate.physicochemical_score
                + 0.35 * (candidate.shape_score or 0)
                + 0.25 * (candidate.esp_score or 0)
            )

        # Sort by combined score
        candidates.sort(key=lambda c: c.combined_score or 0, reverse=True)
        return candidates[: self.top_n_final]

    def get_replacement_molecules(
        self,
        target_fragment: str,
        scaffold_smiles: Optional[str] = None,
    ) -> list[tuple[str, float]]:
        """
        Convenience method to get replacement SMILES with scores.

        Args:
            target_fragment: SMILES of target fragment
            scaffold_smiles: Optional scaffold context

        Returns:
            List of (SMILES, score) tuples
        """
        candidates = self.find_replacements(target_fragment, scaffold_smiles)
        return [(c.smiles, c.combined_score or c.physicochemical_score) for c in candidates]
