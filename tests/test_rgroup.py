"""Tests for RGroup replacement module."""

import pytest
from rdkit import Chem

from expopt.rgroup.replacement import RGroupReplacement, ReplacementCandidate
from expopt.database.mock_database import MockDatabase


class TestRGroupReplacement:
    """Tests for RGroupReplacement class."""

    def test_init_default_database(self):
        """Test initialization with default database."""
        rgroup = RGroupReplacement()
        assert rgroup.database is not None

    def test_init_custom_database(self):
        """Test initialization with custom database."""
        db = MockDatabase()
        rgroup = RGroupReplacement(database=db)
        assert rgroup.database is db

    def test_find_replacements_valid_fragment(self):
        """Test finding replacements for a valid fragment."""
        rgroup = RGroupReplacement()
        target = "NC(=O)C"  # Acetamide

        candidates = rgroup.find_replacements(target)

        assert isinstance(candidates, list)

    def test_find_replacements_with_scaffold(self):
        """Test finding replacements with scaffold context."""
        rgroup = RGroupReplacement()
        target = "NC(=O)C"
        scaffold = "c1ccccc1"

        candidates = rgroup.find_replacements(target, scaffold_smiles=scaffold)

        assert isinstance(candidates, list)

    def test_find_replacements_invalid_smiles(self):
        """Test finding replacements with invalid SMILES."""
        rgroup = RGroupReplacement()

        with pytest.raises(ValueError):
            rgroup.find_replacements("invalid_smiles")

    def test_get_replacement_molecules(self):
        """Test convenience method for getting replacement SMILES."""
        rgroup = RGroupReplacement()
        target = "NC(=O)C"

        results = rgroup.get_replacement_molecules(target)

        assert isinstance(results, list)
        for smiles, score in results:
            assert isinstance(smiles, str)
            assert isinstance(score, float)

    def test_database_search_with_hits(self):
        """Test database search when hits are found."""
        rgroup = RGroupReplacement()
        target = "NC(=O)C"  # This is in the mock database
        target_mol = Chem.MolFromSmiles(target)
        scaffold = "c1ccccc1"

        candidates = rgroup._database_search(target, scaffold, target_mol)

        # Should find alternatives from the mock database
        assert isinstance(candidates, list)

    def test_fallback_generation(self):
        """Test fallback candidate generation."""
        rgroup = RGroupReplacement()
        target_mol = Chem.MolFromSmiles("NC(=O)C")

        candidates = rgroup._fallback_generation(target_mol)

        assert isinstance(candidates, list)

    def test_rank_candidates(self):
        """Test candidate ranking with scores."""
        rgroup = RGroupReplacement()
        target_mol = Chem.MolFromSmiles("NC(=O)C")

        # Create mock candidates
        candidates = [
            ReplacementCandidate(
                smiles="NC(=O)CC",
                mol=Chem.MolFromSmiles("NC(=O)CC"),
                physicochemical_score=0.8,
                source="test",
            ),
            ReplacementCandidate(
                smiles="NC(=O)CCC",
                mol=Chem.MolFromSmiles("NC(=O)CCC"),
                physicochemical_score=0.6,
                source="test",
            ),
        ]

        ranked = rgroup._rank_candidates(candidates, target_mol)

        assert len(ranked) == 2
        # All should have combined scores after ranking
        for candidate in ranked:
            assert candidate.combined_score is not None


class TestReplacementCandidate:
    """Tests for ReplacementCandidate dataclass."""

    def test_create_candidate(self):
        """Test creating a ReplacementCandidate."""
        mol = Chem.MolFromSmiles("CC")
        candidate = ReplacementCandidate(
            smiles="CC",
            mol=mol,
            physicochemical_score=0.9,
            source="database",
        )

        assert candidate.smiles == "CC"
        assert candidate.physicochemical_score == 0.9
        assert candidate.source == "database"
        assert candidate.shape_score is None
        assert candidate.esp_score is None

    def test_candidate_with_all_scores(self):
        """Test candidate with all scores set."""
        mol = Chem.MolFromSmiles("CC")
        candidate = ReplacementCandidate(
            smiles="CC",
            mol=mol,
            physicochemical_score=0.9,
            shape_score=0.8,
            esp_score=0.7,
            combined_score=0.85,
            source="generated",
        )

        assert candidate.shape_score == 0.8
        assert candidate.esp_score == 0.7
        assert candidate.combined_score == 0.85
