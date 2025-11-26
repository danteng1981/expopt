"""
Tests for RGroup Replacement Module
"""

import pytest

try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

from expopt.rgroup import RGroupReplacer, FragmentDatabase, FragmentCandidate


pytestmark = pytest.mark.skipif(not HAS_RDKIT, reason="RDKit not installed")


class TestFragmentDatabase:
    """Tests for FragmentDatabase class."""

    def test_add_fragment(self):
        """Test adding a single fragment."""
        db = FragmentDatabase()
        assert db.add_fragment("c1ccccc1")  # Benzene
        assert len(db) == 1

    def test_add_invalid_fragment(self):
        """Test adding invalid SMILES."""
        db = FragmentDatabase()
        assert not db.add_fragment("invalid_smiles")
        assert len(db) == 0

    def test_add_fragments_from_list(self):
        """Test adding multiple fragments."""
        db = FragmentDatabase()
        smiles_list = ["c1ccccc1", "C1CCCCC1", "c1ccncc1"]  # benzene, cyclohexane, pyridine
        count = db.add_fragments_from_list(smiles_list)
        assert count == 3
        assert len(db) == 3

    def test_lookup_existing(self):
        """Test looking up an existing fragment."""
        db = FragmentDatabase()
        db.add_fragment("c1ccccc1")
        mol = db.lookup("c1ccccc1")
        assert mol is not None

    def test_lookup_missing(self):
        """Test looking up a non-existing fragment."""
        db = FragmentDatabase()
        db.add_fragment("c1ccccc1")
        mol = db.lookup("C1CCCCC1")  # Different molecule
        assert mol is None

    def test_find_similar(self):
        """Test finding similar fragments."""
        db = FragmentDatabase()
        db.add_fragments_from_list([
            "c1ccccc1",      # benzene
            "c1ccc(C)cc1",   # toluene
            "c1ccc(O)cc1",   # phenol
            "C1CCCCC1",      # cyclohexane
        ])

        # Benzene should be most similar to itself and substituted benzenes
        similar = db.find_similar("c1ccccc1", threshold=0.3)
        assert len(similar) > 0
        assert similar[0][1] == 1.0  # Exact match has similarity 1.0


class TestRGroupReplacer:
    """Tests for RGroupReplacer class."""

    def test_exact_match(self):
        """Test exact matching in database."""
        db = FragmentDatabase()
        db.add_fragment("c1ccccc1")

        replacer = RGroupReplacer(database=db)
        result = replacer.exact_match("c1ccccc1")

        assert result is not None
        assert result.source == "database"
        assert result.similarity_score == 1.0

    def test_exact_match_missing(self):
        """Test exact match when fragment not in database."""
        db = FragmentDatabase()
        db.add_fragment("c1ccccc1")

        replacer = RGroupReplacer(database=db)
        result = replacer.exact_match("C1CCCCC1")

        assert result is None

    def test_generate_candidates(self):
        """Test candidate generation."""
        db = FragmentDatabase()
        db.add_fragments_from_list([
            "c1ccccc1",      # benzene
            "c1ccc(C)cc1",   # toluene
            "c1ccc(O)cc1",   # phenol
            "c1ccc(N)cc1",   # aniline
            "C1CCCCC1",      # cyclohexane (less similar)
        ])

        replacer = RGroupReplacer(database=db)
        candidates = replacer.generate_candidates(
            "c1ccc(F)cc1",  # fluorobenzene
            similarity_threshold=0.3,  # Lower threshold for topological fingerprints
            max_candidates=5
        )

        assert len(candidates) > 0
        # All candidates should be reasonably similar
        for c in candidates:
            assert c.similarity_score >= 0.3

    def test_find_replacements_priority(self):
        """Test that exact match is prioritized."""
        db = FragmentDatabase()
        db.add_fragments_from_list([
            "c1ccccc1",      # benzene (exact match)
            "c1ccc(C)cc1",   # toluene
        ])

        replacer = RGroupReplacer(database=db)
        results = replacer.find_replacements(
            "c1ccccc1",
            similarity_threshold=0.3
        )

        assert len(results) > 0
        # First result should be exact match
        assert results[0].source == "database"
        assert results[0].similarity_score == 1.0

    def test_property_constraints(self):
        """Test filtering by property constraints."""
        db = FragmentDatabase()
        db.add_fragments_from_list([
            "c1ccccc1",           # MW ~78
            "c1ccc(CCCCC)cc1",    # MW ~148 (larger)
            "c1ccc(CC)cc1",       # MW ~106
        ])

        replacer = RGroupReplacer(database=db)
        candidates = replacer.generate_candidates(
            "c1ccccc1",
            similarity_threshold=0.3,
            property_constraints={"mw": (70, 120)}  # Limit MW
        )

        # Should filter out the large molecule
        for c in candidates:
            assert c.properties["mw"] <= 120

    def test_get_attachment_points(self):
        """Test finding attachment points."""
        replacer = RGroupReplacer()

        # Fragment with dummy atom attachment point
        points = replacer.get_attachment_points("*c1ccccc1")
        # Should find the attachment point
        assert len(points) > 0
