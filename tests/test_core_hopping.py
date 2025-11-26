"""
Tests for Core Hopping Module
"""

import pytest

try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

from expopt.core_hopping import CoreHopper, CoreDatabase, CoreCandidate


pytestmark = pytest.mark.skipif(not HAS_RDKIT, reason="RDKit not installed")


class TestCoreDatabase:
    """Tests for CoreDatabase class."""

    def test_add_core(self):
        """Test adding a single core."""
        db = CoreDatabase()
        assert db.add_core("c1ccc2ccccc2c1")  # naphthalene
        assert len(db) == 1

    def test_add_invalid_core(self):
        """Test adding invalid SMILES."""
        db = CoreDatabase()
        assert not db.add_core("invalid_smiles")
        assert len(db) == 0

    def test_add_cores_from_list(self):
        """Test adding multiple cores."""
        db = CoreDatabase()
        smiles_list = [
            "c1ccc2ccccc2c1",  # naphthalene
            "c1ccc2ncccc2c1",  # quinoline
            "c1ccc2[nH]ccc2c1",  # indole
        ]
        count = db.add_cores_from_list(smiles_list)
        assert count == 3
        assert len(db) == 3

    def test_lookup_existing(self):
        """Test looking up an existing core."""
        db = CoreDatabase()
        db.add_core("c1ccc2ccccc2c1")
        mol = db.lookup("c1ccc2ccccc2c1")
        assert mol is not None

    def test_find_similar_cores(self):
        """Test finding similar cores."""
        db = CoreDatabase()
        db.add_cores_from_list([
            "c1ccc2ccccc2c1",   # naphthalene
            "c1ccc2ncccc2c1",   # quinoline
            "c1ccc2[nH]ccc2c1", # indole
            "c1ccccc1",         # benzene (smaller, less similar)
        ])

        similar = db.find_similar("c1ccc2ccccc2c1", threshold=0.3)
        assert len(similar) > 0
        assert similar[0][1] == 1.0  # Exact match


class TestCoreHopper:
    """Tests for CoreHopper class."""

    def test_extract_core(self):
        """Test extracting Murcko scaffold."""
        hopper = CoreHopper()

        # Aspirin-like structure
        core = hopper.extract_core("CC(=O)Oc1ccccc1C(=O)O")
        assert core is not None
        # Should extract the benzene ring scaffold
        assert "c1ccccc1" in core or Chem.MolFromSmiles(core).GetRingInfo().NumRings() >= 1

    def test_extract_core_simple(self):
        """Test core extraction on simple molecule."""
        hopper = CoreHopper()
        core = hopper.extract_core("c1ccc(C)cc1")  # toluene
        assert core is not None
        # Core should be benzene
        assert Chem.MolFromSmiles(core).GetRingInfo().NumRings() == 1

    def test_extract_generic_core(self):
        """Test generic scaffold extraction."""
        hopper = CoreHopper()
        generic = hopper.extract_generic_core("c1ncc(N)nc1")  # pyrimidine with amino
        assert generic is not None
        # Generic core should have only carbons
        mol = Chem.MolFromSmiles(generic)
        assert mol is not None

    def test_decompose_molecule(self):
        """Test molecule decomposition."""
        hopper = CoreHopper()
        result = hopper.decompose_molecule("c1ccc(CC)cc1")  # ethylbenzene

        assert result is not None
        assert "core" in result
        assert "core_atoms" in result
        assert "rgroup_atoms" in result

    def test_find_hopping_candidates(self):
        """Test finding scaffold hopping candidates."""
        db = CoreDatabase()
        db.add_cores_from_list([
            "c1ccc2ccccc2c1",   # naphthalene
            "c1ccc2ncccc2c1",   # quinoline
            "c1ccc2[nH]ccc2c1", # indole
            "c1ccc2occc2c1",    # benzofuran
        ])

        hopper = CoreHopper(database=db)
        candidates = hopper.find_hopping_candidates(
            "c1ccc2ccccc2c1C",  # methylnaphthalene
            similarity_threshold=0.3,
            max_candidates=10
        )

        assert len(candidates) > 0
        # Should find naphthalene and similar bicyclic structures
        for c in candidates:
            assert c.similarity_score >= 0.3

    def test_get_compatible_cores(self):
        """Test finding cores with required attachment points."""
        db = CoreDatabase()
        db.add_cores_from_list([
            "c1ccccc1",         # benzene - 6 potential attachments
            "c1ccc2ccccc2c1",   # naphthalene - 8 potential attachments
        ])

        hopper = CoreHopper(database=db)
        compatible = hopper.get_compatible_cores(
            "c1ccccc1C",  # toluene
            required_attachment_count=3
        )

        # All compatible cores should have >= 3 attachment points
        for c in compatible:
            assert len(c.attachment_points) >= 3

    def test_calculate_hopping_potential(self):
        """Test hopping potential calculation."""
        hopper = CoreHopper()
        potential = hopper.calculate_hopping_potential(
            "c1ccc2ccccc2c1",  # naphthalene
            "c1ccc2ncccc2c1"   # quinoline
        )

        assert "similarity" in potential
        assert "novelty" in potential
        assert "mw_diff" in potential
        assert potential["novelty"] == 1.0 - potential["similarity"]


class TestCoreHopperEdgeCases:
    """Edge case tests for CoreHopper."""

    def test_extract_core_invalid_smiles(self):
        """Test core extraction with invalid SMILES."""
        hopper = CoreHopper()
        core = hopper.extract_core("invalid")
        assert core is None

    def test_decompose_invalid_molecule(self):
        """Test decomposition with invalid SMILES."""
        hopper = CoreHopper()
        result = hopper.decompose_molecule("invalid")
        assert result is None

    def test_empty_database(self):
        """Test with empty database."""
        hopper = CoreHopper()
        candidates = hopper.find_hopping_candidates("c1ccccc1C")
        assert candidates == []
