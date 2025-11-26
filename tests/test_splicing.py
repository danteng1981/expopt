"""
Tests for Splicing Module
"""

import pytest

try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

from expopt.splicing import RGroupSplicer, CoreSplicer, SplicedMolecule


pytestmark = pytest.mark.skipif(not HAS_RDKIT, reason="RDKit not installed")


class TestRGroupSplicer:
    """Tests for RGroupSplicer class."""

    def test_basic_splice(self):
        """Test basic splicing of two molecules."""
        splicer = RGroupSplicer()

        # Splice methyl onto benzene
        result = splicer.splice(
            "c1ccccc1",  # benzene
            "C",         # methyl
            0,           # attach to first carbon of benzene
            0            # attach to carbon of methyl
        )

        assert result is not None
        assert result.smiles is not None
        # Result should be toluene-like
        mol = Chem.MolFromSmiles(result.smiles)
        assert mol is not None
        assert mol.GetNumHeavyAtoms() == 7  # 6 from benzene + 1 from methyl

    def test_splice_with_dummy_atoms(self):
        """Test splicing using dummy atoms as attachment points."""
        splicer = RGroupSplicer()

        # Molecules with dummy atoms marking attachment
        result = splicer.splice_with_dummy_atoms(
            "*c1ccccc1",  # benzene with attachment point
            "*C"          # methyl with attachment point
        )

        assert result is not None
        mol = Chem.MolFromSmiles(result.smiles)
        assert mol is not None
        # Should have no dummy atoms in result
        for atom in mol.GetAtoms():
            assert atom.GetAtomicNum() != 0

    def test_splice_multiple_fragments(self):
        """Test splicing multiple fragments."""
        splicer = RGroupSplicer()

        # Splice two methyls onto benzene at different positions
        result = splicer.splice_multiple_fragments(
            "c1ccccc1",  # benzene
            [
                ("C", 0, 0),  # methyl at position 0
                ("C", 3, 0),  # methyl at position 3 (para)
            ]
        )

        assert result is not None
        mol = Chem.MolFromSmiles(result.smiles)
        assert mol is not None
        assert mol.GetNumHeavyAtoms() == 8  # 6 + 2

    def test_splice_invalid_molecules(self):
        """Test splicing with invalid molecules."""
        splicer = RGroupSplicer()

        result = splicer.splice(
            "invalid",
            "C",
            0, 0
        )
        assert result is None

        result = splicer.splice(
            "c1ccccc1",
            "invalid",
            0, 0
        )
        assert result is None


class TestCoreSplicer:
    """Tests for CoreSplicer class."""

    def test_enumerate_connections(self):
        """Test enumeration of connection combinations."""
        splicer = CoreSplicer()

        # Simple case: benzene with one R-group
        connections = list(splicer.enumerate_connections(
            "c1ccccc1",  # benzene
            ["C"],       # methyl
            core_attachment_points=[0, 1, 2]  # limit attachment points
        ))

        assert len(connections) > 0
        # Should enumerate different attachment possibilities

    def test_splice_core_with_rgroups(self):
        """Test core splicing with R-groups."""
        splicer = CoreSplicer()

        result = splicer.splice_core_with_rgroups(
            "c1ccccc1",   # benzene core
            ["C", "O"],   # methyl and hydroxyl
            [(0, 0), (3, 0)]  # connections at positions 0 and 3
        )

        assert result is not None
        assert result.smiles is not None
        mol = Chem.MolFromSmiles(result.smiles)
        assert mol is not None

    def test_find_best_connections(self):
        """Test finding best connections by similarity."""
        splicer = CoreSplicer()

        # Find best way to attach methyls to benzene to be similar to xylene
        results = splicer.find_best_connections(
            "c1ccccc1",       # benzene core
            ["C", "C"],       # two methyls
            "Cc1ccc(C)cc1",   # p-xylene target
            max_results=5,
            core_attachment_points=[0, 1, 2, 3, 4, 5]
        )

        assert len(results) > 0
        # Best result should have reasonable similarity
        assert results[0].similarity_to_target > 0

    def test_score_molecule(self):
        """Test molecule scoring by similarity."""
        splicer = CoreSplicer()

        # Identical molecules should have similarity 1.0
        score = splicer.score_molecule("c1ccccc1", "c1ccccc1")
        assert score == 1.0

        # Similar molecules should have some similarity (not 0)
        score = splicer.score_molecule("c1ccc(C)cc1", "c1ccccc1")
        assert score > 0.1  # Toluene and benzene should have some similarity

        # Invalid SMILES should return 0
        score = splicer.score_molecule("invalid", "c1ccccc1")
        assert score == 0.0

    def test_batch_splice_and_score(self):
        """Test batch processing of multiple cores."""
        splicer = CoreSplicer()

        cores = [
            "c1ccccc1",      # benzene
            "c1ccncc1",      # pyridine
        ]
        rgroups = ["C"]     # methyl
        target = "c1ccc(C)cc1"  # toluene

        results = splicer.batch_splice_and_score(
            cores, rgroups, target, top_n=5
        )

        assert len(results) > 0
        # Results should be sorted by similarity
        for i in range(len(results) - 1):
            assert results[i].similarity_to_target >= results[i + 1].similarity_to_target


class TestCoreSplicerEdgeCases:
    """Edge case tests for CoreSplicer."""

    def test_empty_rgroups(self):
        """Test with empty R-groups list."""
        splicer = CoreSplicer()

        results = splicer.find_best_connections(
            "c1ccccc1",
            [],
            "c1ccccc1",
            max_results=5
        )
        # Should handle gracefully
        assert isinstance(results, list)

    def test_mismatched_connections(self):
        """Test with mismatched connections and R-groups."""
        splicer = CoreSplicer()

        result = splicer.splice_core_with_rgroups(
            "c1ccccc1",
            ["C"],           # 1 R-group
            [(0, 0), (1, 0)] # 2 connections - mismatch
        )
        assert result is None

    def test_invalid_core(self):
        """Test with invalid core SMILES."""
        splicer = CoreSplicer()

        result = splicer.splice_core_with_rgroups(
            "invalid",
            ["C"],
            [(0, 0)]
        )
        assert result is None


class TestSplicedMolecule:
    """Tests for SplicedMolecule dataclass."""

    def test_creation(self):
        """Test creating SplicedMolecule."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        spliced = SplicedMolecule(
            smiles="c1ccccc1",
            mol=mol,
            similarity_to_target=0.5
        )

        assert spliced.smiles == "c1ccccc1"
        assert spliced.similarity_to_target == 0.5
        assert spliced.connection_info == {}

    def test_default_values(self):
        """Test default values in SplicedMolecule."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        spliced = SplicedMolecule(smiles="c1ccccc1", mol=mol)

        assert spliced.similarity_to_target == 0.0
        assert spliced.connection_info == {}
