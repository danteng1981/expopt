"""Tests for CoreHopping module."""

import pytest
from rdkit import Chem

from expopt.core_hopping.hopping import CoreHopping, CoreCandidate
from expopt.database.mock_database import MockDatabase


class TestCoreHopping:
    """Tests for CoreHopping class."""

    def test_init_default_database(self):
        """Test initialization with default database."""
        hopping = CoreHopping()
        assert hopping.database is not None

    def test_init_custom_database(self):
        """Test initialization with custom database."""
        db = MockDatabase()
        hopping = CoreHopping(database=db)
        assert hopping.database is db

    def test_find_alternative_cores_benzene(self):
        """Test finding alternative cores for benzene."""
        hopping = CoreHopping()
        target = "c1ccccc1"  # Benzene

        candidates = hopping.find_alternative_cores(target)

        assert isinstance(candidates, list)

    def test_find_alternative_cores_with_min_sites(self):
        """Test finding cores with minimum connection sites."""
        hopping = CoreHopping()
        target = "c1ccccc1"

        candidates = hopping.find_alternative_cores(target, min_connection_sites=4)

        assert isinstance(candidates, list)
        for candidate in candidates:
            assert candidate.connection_sites >= 4

    def test_find_alternative_cores_invalid_smiles(self):
        """Test finding alternatives with invalid SMILES."""
        hopping = CoreHopping()

        with pytest.raises(ValueError):
            hopping.find_alternative_cores("invalid_smiles")

    def test_get_alternative_core_molecules(self):
        """Test convenience method for getting alternative core SMILES."""
        hopping = CoreHopping()
        target = "c1ccccc1"

        results = hopping.get_alternative_core_molecules(target)

        assert isinstance(results, list)
        for smiles, score in results:
            assert isinstance(smiles, str)
            assert isinstance(score, float)

    def test_database_search(self):
        """Test database search for cores."""
        hopping = CoreHopping()
        target = "c1ccccc1"
        target_mol = Chem.MolFromSmiles(target)

        candidates = hopping._database_search(target, target_mol, min_sites=1)

        assert isinstance(candidates, list)

    def test_fallback_generation(self):
        """Test fallback core generation."""
        hopping = CoreHopping()
        target_mol = Chem.MolFromSmiles("c1ccccc1")

        candidates = hopping._fallback_generation(target_mol, min_sites=4)

        assert isinstance(candidates, list)

    def test_apply_constraints(self):
        """Test constraint application."""
        hopping = CoreHopping()
        target_mol = Chem.MolFromSmiles("c1ccccc1")

        # Create mock candidates
        candidates = [
            CoreCandidate(
                smiles="c1ccncc1",
                mol=Chem.MolFromSmiles("c1ccncc1"),
                connection_sites=5,
                physicochemical_score=0.8,
            ),
            CoreCandidate(
                smiles="c1ccc2ccccc2c1",
                mol=Chem.MolFromSmiles("c1ccc2ccccc2c1"),
                connection_sites=8,
                physicochemical_score=0.6,
            ),
        ]

        filtered = hopping._apply_constraints(candidates, target_mol, min_sites=5)

        assert all(c.connection_sites >= 5 for c in filtered)

    def test_rank_candidates(self):
        """Test candidate ranking."""
        hopping = CoreHopping()
        target_mol = Chem.MolFromSmiles("c1ccccc1")

        candidates = [
            CoreCandidate(
                smiles="c1ccncc1",
                mol=Chem.MolFromSmiles("c1ccncc1"),
                connection_sites=5,
                physicochemical_score=0.8,
            ),
            CoreCandidate(
                smiles="c1ccoc1",
                mol=Chem.MolFromSmiles("c1ccoc1"),
                connection_sites=4,
                physicochemical_score=0.7,
            ),
        ]

        ranked = hopping._rank_candidates(candidates, target_mol)

        assert len(ranked) == 2
        for candidate in ranked:
            assert candidate.combined_score is not None


class TestCoreCandidate:
    """Tests for CoreCandidate dataclass."""

    def test_create_candidate(self):
        """Test creating a CoreCandidate."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        candidate = CoreCandidate(
            smiles="c1ccccc1",
            mol=mol,
            connection_sites=6,
            physicochemical_score=0.9,
            source="database",
        )

        assert candidate.smiles == "c1ccccc1"
        assert candidate.connection_sites == 6
        assert candidate.physicochemical_score == 0.9
        assert candidate.source == "database"

    def test_candidate_with_all_scores(self):
        """Test candidate with all scores."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        candidate = CoreCandidate(
            smiles="c1ccccc1",
            mol=mol,
            connection_sites=6,
            physicochemical_score=0.9,
            shape_score=0.85,
            esp_score=0.8,
            combined_score=0.87,
            source="generated",
        )

        assert candidate.shape_score == 0.85
        assert candidate.esp_score == 0.8
        assert candidate.combined_score == 0.87
