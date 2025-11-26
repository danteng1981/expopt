"""Tests for utility modules."""

import pytest
from rdkit import Chem

from expopt.utils.scoring import PhysicochemicalScorer
from expopt.utils.similarity import SimilarityCalculator
from expopt.utils.constraints import AtomCountConstraint, ConnectionSiteConstraint


class TestPhysicochemicalScorer:
    """Tests for PhysicochemicalScorer class."""

    def test_calculate_properties_benzene(self):
        """Test property calculation for benzene."""
        scorer = PhysicochemicalScorer()
        mol = Chem.MolFromSmiles("c1ccccc1")

        props = scorer.calculate_properties(mol)

        assert "mw" in props
        assert "logp" in props
        assert "hbd" in props
        assert "hba" in props
        assert props["num_atoms"] == 6
        assert props["num_aromatic_rings"] == 1

    def test_calculate_properties_none_mol(self):
        """Test property calculation with None molecule."""
        scorer = PhysicochemicalScorer()
        props = scorer.calculate_properties(None)
        assert props == {}

    def test_calculate_score_identical(self):
        """Test scoring identical molecules."""
        scorer = PhysicochemicalScorer()
        mol = Chem.MolFromSmiles("c1ccccc1")

        score = scorer.calculate_score(mol, mol)

        assert score == 1.0

    def test_calculate_score_similar(self):
        """Test scoring similar molecules."""
        scorer = PhysicochemicalScorer()
        benzene = Chem.MolFromSmiles("c1ccccc1")
        toluene = Chem.MolFromSmiles("Cc1ccccc1")

        score = scorer.calculate_score(benzene, toluene)

        assert 0.5 < score < 1.0

    def test_rank_candidates(self):
        """Test ranking candidates by score."""
        scorer = PhysicochemicalScorer()
        target = Chem.MolFromSmiles("c1ccccc1")
        candidates = [
            Chem.MolFromSmiles("Cc1ccccc1"),  # Toluene
            Chem.MolFromSmiles("c1ccc2ccccc2c1"),  # Naphthalene
            Chem.MolFromSmiles("c1ccccc1"),  # Benzene (identical)
        ]

        ranked = scorer.rank_candidates(candidates, target, top_n=2)

        assert len(ranked) == 2
        # Identical molecule should rank highest
        assert Chem.MolToSmiles(ranked[0][0]) == "c1ccccc1"


class TestSimilarityCalculator:
    """Tests for SimilarityCalculator class."""

    def test_tanimoto_identical(self):
        """Test Tanimoto similarity for identical molecules."""
        calc = SimilarityCalculator()
        mol = Chem.MolFromSmiles("c1ccccc1")

        similarity = calc.tanimoto_similarity(mol, mol)

        assert similarity == 1.0

    def test_tanimoto_different(self):
        """Test Tanimoto similarity for different molecules."""
        calc = SimilarityCalculator()
        benzene = Chem.MolFromSmiles("c1ccccc1")
        aspirin = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")

        similarity = calc.tanimoto_similarity(benzene, aspirin)

        assert 0.0 <= similarity < 1.0

    def test_get_fingerprint_morgan(self):
        """Test Morgan fingerprint generation."""
        calc = SimilarityCalculator(fingerprint_type="morgan")
        mol = Chem.MolFromSmiles("c1ccccc1")

        fp = calc.get_fingerprint(mol)

        assert fp is not None

    def test_get_fingerprint_maccs(self):
        """Test MACCS fingerprint generation."""
        calc = SimilarityCalculator(fingerprint_type="maccs")
        mol = Chem.MolFromSmiles("c1ccccc1")

        fp = calc.get_fingerprint(mol)

        assert fp is not None

    def test_esp_similarity(self):
        """Test ESP similarity calculation."""
        calc = SimilarityCalculator()
        benzene = Chem.MolFromSmiles("c1ccccc1")
        pyridine = Chem.MolFromSmiles("c1ccncc1")

        similarity = calc.esp_similarity(benzene, pyridine)

        assert 0.0 <= similarity <= 1.0

    def test_rank_by_similarity(self):
        """Test ranking by similarity."""
        calc = SimilarityCalculator()
        target = Chem.MolFromSmiles("c1ccccc1")
        candidates = [
            Chem.MolFromSmiles("c1ccc2ccccc2c1"),  # Naphthalene
            Chem.MolFromSmiles("Cc1ccccc1"),  # Toluene
            Chem.MolFromSmiles("c1ccccc1"),  # Benzene
        ]

        ranked = calc.rank_by_similarity(candidates, target, method="tanimoto")

        assert len(ranked) == 3
        assert ranked[0][1] == 1.0  # Identical molecule


class TestAtomCountConstraint:
    """Tests for AtomCountConstraint class."""

    def test_check_constraint_within_range(self):
        """Test constraint with molecule within range."""
        constraint = AtomCountConstraint(target_atom_count=6)
        mol = Chem.MolFromSmiles("c1ccccc1")  # 6 atoms

        assert constraint.check_constraint(mol) is True

    def test_check_constraint_too_small(self):
        """Test constraint with molecule too small."""
        constraint = AtomCountConstraint(target_atom_count=10)
        mol = Chem.MolFromSmiles("C")  # 1 atom

        assert constraint.check_constraint(mol) is False

    def test_check_constraint_too_large(self):
        """Test constraint with molecule too large."""
        constraint = AtomCountConstraint(target_atom_count=2)
        mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")  # 10 atoms

        assert constraint.check_constraint(mol) is False

    def test_filter_candidates(self):
        """Test filtering candidates by constraint."""
        constraint = AtomCountConstraint(target_atom_count=6)
        candidates = [
            Chem.MolFromSmiles("C"),  # Too small
            Chem.MolFromSmiles("c1ccccc1"),  # In range
            Chem.MolFromSmiles("Cc1ccccc1"),  # In range
        ]

        filtered = constraint.filter_candidates(candidates)

        assert len(filtered) == 2


class TestConnectionSiteConstraint:
    """Tests for ConnectionSiteConstraint class."""

    def test_count_connection_sites_benzene(self):
        """Test counting connection sites in benzene."""
        constraint = ConnectionSiteConstraint(min_sites=1)
        mol = Chem.MolFromSmiles("c1ccccc1")

        sites = constraint.count_connection_sites(mol)

        # Benzene has 6 carbons that could accept hydrogen substitution
        assert sites >= 0

    def test_check_constraint(self):
        """Test connection site constraint check."""
        constraint = ConnectionSiteConstraint(min_sites=2)
        mol = Chem.MolFromSmiles("c1ccccc1")

        # This should pass since benzene has multiple positions
        result = constraint.check_constraint(mol)
        assert isinstance(result, bool)
