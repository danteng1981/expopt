"""Tests for database module."""

import pytest
from rdkit import Chem

from expopt.database.interface import (
    DatabaseInterface,
    PatentRecord,
    ScaffoldRecord,
    CoreRecord,
)
from expopt.database.mock_database import MockDatabase


class TestPatentRecord:
    """Tests for PatentRecord dataclass."""

    def test_create_patent_record(self):
        """Test creating a PatentRecord."""
        record = PatentRecord(
            patent_id="PAT001",
            smiles="c1ccccc1",
            scaffold_smiles="c1ccccc1",
            fragments=["C", "CC"],
        )

        assert record.patent_id == "PAT001"
        assert record.smiles == "c1ccccc1"
        assert record.fragments == ["C", "CC"]


class TestScaffoldRecord:
    """Tests for ScaffoldRecord dataclass."""

    def test_create_scaffold_record(self):
        """Test creating a ScaffoldRecord."""
        record = ScaffoldRecord(
            scaffold_id="SCF001",
            scaffold_smiles="c1ccccc1",
            fragments=["C", "CC"],
            patent_ids=["PAT001"],
        )

        assert record.scaffold_id == "SCF001"
        assert record.scaffold_smiles == "c1ccccc1"


class TestCoreRecord:
    """Tests for CoreRecord dataclass."""

    def test_create_core_record(self):
        """Test creating a CoreRecord."""
        record = CoreRecord(
            core_id="CORE001",
            core_smiles="c1ccccc1",
            connection_sites=6,
            patent_ids=["PAT001"],
        )

        assert record.core_id == "CORE001"
        assert record.connection_sites == 6


class TestMockDatabase:
    """Tests for MockDatabase class."""

    def test_init(self):
        """Test MockDatabase initialization."""
        db = MockDatabase()
        assert db is not None

    def test_search_by_fragment(self):
        """Test searching by fragment."""
        db = MockDatabase()

        results = db.search_by_fragment("NC(=O)C")

        assert isinstance(results, list)

    def test_search_by_scaffold(self):
        """Test searching by scaffold."""
        db = MockDatabase()

        results = db.search_by_scaffold("c1ccccc1")

        assert isinstance(results, list)
        assert len(results) > 0

    def test_search_by_scaffold_not_found(self):
        """Test searching for non-existent scaffold."""
        db = MockDatabase()

        results = db.search_by_scaffold("nonexistent")

        assert results == []

    def test_get_alternative_fragments(self):
        """Test getting alternative fragments."""
        db = MockDatabase()

        fragments = db.get_alternative_fragments("c1ccccc1")

        assert isinstance(fragments, list)
        assert len(fragments) > 0

    def test_get_alternative_fragments_not_found(self):
        """Test getting alternatives for non-existent scaffold."""
        db = MockDatabase()

        fragments = db.get_alternative_fragments("nonexistent")

        assert fragments == []

    def test_search_cores(self):
        """Test searching for cores."""
        db = MockDatabase()

        cores = db.search_cores(min_sites=1)

        assert isinstance(cores, list)
        assert len(cores) > 0
        assert all(c.connection_sites >= 1 for c in cores)

    def test_search_cores_min_sites(self):
        """Test searching cores with minimum sites."""
        db = MockDatabase()

        cores = db.search_cores(min_sites=6)

        assert all(c.connection_sites >= 6 for c in cores)

    def test_get_alternative_cores(self):
        """Test getting alternative cores."""
        db = MockDatabase()

        alternatives = db.get_alternative_cores("c1ccccc1")

        assert isinstance(alternatives, list)
        # Should not include the target core itself
        for core in alternatives:
            assert core.core_smiles != "c1ccccc1"

    def test_get_fragment_library(self):
        """Test getting fragment library."""
        db = MockDatabase()

        library = db.get_fragment_library()

        assert isinstance(library, list)
        assert len(library) > 0

    def test_add_fragment(self):
        """Test adding a fragment."""
        db = MockDatabase()
        initial_count = len(db.get_fragment_library())

        db.add_fragment("CCCC")

        assert len(db.get_fragment_library()) == initial_count + 1

    def test_add_fragment_duplicate(self):
        """Test adding duplicate fragment."""
        db = MockDatabase()
        db.add_fragment("CCCC")
        count_after_first = len(db.get_fragment_library())

        db.add_fragment("CCCC")

        assert len(db.get_fragment_library()) == count_after_first

    def test_add_core(self):
        """Test adding a core."""
        db = MockDatabase()
        initial_count = len(db.search_cores())

        db.add_core("c1ccc2ncccc2c1", 7, ["PAT999"])

        assert len(db.search_cores()) == initial_count + 1

    def test_smiles_to_mol(self):
        """Test SMILES to Mol conversion."""
        db = MockDatabase()

        mol = db.smiles_to_mol("c1ccccc1")

        assert mol is not None
        assert mol.GetNumAtoms() == 6

    def test_smiles_to_mol_invalid(self):
        """Test SMILES to Mol with invalid SMILES."""
        db = MockDatabase()

        mol = db.smiles_to_mol("invalid")

        assert mol is None

    def test_mol_to_smiles(self):
        """Test Mol to SMILES conversion."""
        db = MockDatabase()
        mol = Chem.MolFromSmiles("c1ccccc1")

        smiles = db.mol_to_smiles(mol)

        assert smiles is not None
        assert isinstance(smiles, str)

    def test_mol_to_smiles_none(self):
        """Test Mol to SMILES with None input."""
        db = MockDatabase()

        smiles = db.mol_to_smiles(None)

        assert smiles is None
