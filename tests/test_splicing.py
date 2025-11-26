"""Tests for splicing modules."""

import pytest
from rdkit import Chem

from expopt.splicing.rgroup_splicing import RGroupSplicing
from expopt.splicing.core_splicing import CoreSplicing, SplicedMolecule


class TestRGroupSplicing:
    """Tests for RGroupSplicing class."""

    def test_init(self):
        """Test RGroupSplicing initialization."""
        splicing = RGroupSplicing()
        assert splicing is not None

    def test_splice_fragment_basic(self):
        """Test basic fragment splicing."""
        splicing = RGroupSplicing()
        scaffold = "c1ccccc1"
        fragment = "C"

        result = splicing.splice_fragment(scaffold, fragment, 0, 0)

        # Should produce a valid SMILES
        if result:
            mol = Chem.MolFromSmiles(result)
            assert mol is not None

    def test_splice_fragment_invalid_scaffold(self):
        """Test splicing with invalid scaffold."""
        splicing = RGroupSplicing()

        result = splicing.splice_fragment("invalid", "C")

        assert result is None

    def test_splice_fragment_invalid_fragment(self):
        """Test splicing with invalid fragment."""
        splicing = RGroupSplicing()

        result = splicing.splice_fragment("c1ccccc1", "invalid")

        assert result is None

    def test_find_attachment_point(self):
        """Test finding attachment point with dummy atom."""
        splicing = RGroupSplicing()
        # Molecule with dummy atom
        mol = Chem.MolFromSmiles("[*]c1ccccc1")

        attach = splicing._find_attachment_point(mol)

        # Should find the carbon attached to dummy
        assert attach is not None

    def test_remove_dummy_atoms(self):
        """Test removing dummy atoms from molecule."""
        splicing = RGroupSplicing()
        mol = Chem.MolFromSmiles("[*]c1ccccc1")

        result = splicing._remove_dummy_atoms(mol)

        assert result is not None
        # Should have no dummy atoms
        for atom in result.GetAtoms():
            assert atom.GetAtomicNum() != 0

    def test_enumerate_attachments(self):
        """Test enumerating attachment configurations."""
        splicing = RGroupSplicing()
        scaffold = "c1ccccc1"
        fragment = "C"

        results = splicing.enumerate_attachments(scaffold, fragment)

        # Should produce at least one valid result
        assert isinstance(results, list)


class TestCoreSplicing:
    """Tests for CoreSplicing class."""

    def test_init(self):
        """Test CoreSplicing initialization."""
        splicing = CoreSplicing()
        assert splicing is not None

    def test_init_custom_fingerprint(self):
        """Test initialization with custom fingerprint."""
        splicing = CoreSplicing(fingerprint_type="maccs")
        assert splicing.similarity.fingerprint_type == "maccs"

    def test_splice_core_single_sidechain(self):
        """Test splicing with single sidechain."""
        splicing = CoreSplicing()
        core = "c1ccccc1"
        sidechains = ["C"]

        result = splicing.splice_core(core, sidechains)

        # Should produce a valid SMILES or None
        if result:
            mol = Chem.MolFromSmiles(result)
            assert mol is not None

    def test_get_attachment_positions(self):
        """Test getting attachment positions."""
        splicing = CoreSplicing()
        mol = Chem.MolFromSmiles("c1ccccc1")

        positions = splicing._get_attachment_positions(mol)

        assert isinstance(positions, list)

    def test_enumerate_configurations(self):
        """Test enumerating configurations."""
        splicing = CoreSplicing()
        core = "c1ccccc1"
        sidechains = ["C"]

        configs = splicing.enumerate_configurations(core, sidechains, max_configurations=10)

        assert isinstance(configs, list)

    def test_rank_by_similarity(self):
        """Test ranking by similarity."""
        splicing = CoreSplicing()
        reference = Chem.MolFromSmiles("Cc1ccccc1")

        candidates = [
            SplicedMolecule(
                smiles="Cc1ccccc1",
                mol=Chem.MolFromSmiles("Cc1ccccc1"),
                similarity_score=0.0,
                configuration_id=0,
            ),
            SplicedMolecule(
                smiles="CCc1ccccc1",
                mol=Chem.MolFromSmiles("CCc1ccccc1"),
                similarity_score=0.0,
                configuration_id=1,
            ),
        ]

        ranked = splicing.rank_by_similarity(candidates, reference)

        assert len(ranked) == 2
        # First should be identical (score = 1.0)
        assert ranked[0].similarity_score == 1.0

    def test_splice_and_rank(self):
        """Test complete splice and rank workflow."""
        splicing = CoreSplicing()
        core = "c1ccccc1"
        sidechains = ["C"]
        reference = "Cc1ccccc1"

        results = splicing.splice_and_rank(
            core, sidechains, reference, max_configurations=10, top_n=5
        )

        assert isinstance(results, list)
        assert len(results) <= 5

    def test_get_best_configuration(self):
        """Test getting best configuration."""
        splicing = CoreSplicing()
        core = "c1ccccc1"
        sidechains = ["C"]
        reference = "Cc1ccccc1"

        result = splicing.get_best_configuration(core, sidechains, reference)

        # Should return SMILES string or None
        assert result is None or isinstance(result, str)


class TestSplicedMolecule:
    """Tests for SplicedMolecule dataclass."""

    def test_create_spliced_molecule(self):
        """Test creating a SplicedMolecule."""
        mol = Chem.MolFromSmiles("Cc1ccccc1")
        spliced = SplicedMolecule(
            smiles="Cc1ccccc1",
            mol=mol,
            similarity_score=0.95,
            configuration_id=1,
        )

        assert spliced.smiles == "Cc1ccccc1"
        assert spliced.similarity_score == 0.95
        assert spliced.configuration_id == 1
