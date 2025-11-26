"""Mock database implementation for testing and demonstration."""

from typing import Optional

from expopt.database.interface import (
    CoreRecord,
    DatabaseInterface,
    PatentRecord,
    ScaffoldRecord,
)


class MockDatabase(DatabaseInterface):
    """
    Mock database implementation for testing purposes.

    This provides sample data for demonstrating the RGroup and CoreHopping workflows.
    In production, this would be replaced with a real database connection.
    """

    def __init__(self):
        """Initialize mock database with sample data."""
        # Sample patent records
        self._patents = {
            "PAT001": PatentRecord(
                patent_id="PAT001",
                smiles="c1ccc(NC(=O)c2ccccc2)cc1",
                scaffold_smiles="c1ccccc1",
                fragments=["NC(=O)c1ccccc1", "NC(=O)C"],
            ),
            "PAT002": PatentRecord(
                patent_id="PAT002",
                smiles="c1ccc(NC(=O)C)cc1",
                scaffold_smiles="c1ccccc1",
                fragments=["NC(=O)C", "NC(=O)CC"],
            ),
            "PAT003": PatentRecord(
                patent_id="PAT003",
                smiles="c1ccc2c(c1)CCN2",
                scaffold_smiles="c1ccc2c(c1)CCN2",
                fragments=["C", "CC"],
            ),
        }

        # Sample scaffold records
        self._scaffolds = {
            "c1ccccc1": ScaffoldRecord(
                scaffold_id="SCF001",
                scaffold_smiles="c1ccccc1",
                fragments=["NC(=O)c1ccccc1", "NC(=O)C", "NC(=O)CC", "F", "Cl", "Br"],
                patent_ids=["PAT001", "PAT002"],
            ),
            "c1ccc2c(c1)CCN2": ScaffoldRecord(
                scaffold_id="SCF002",
                scaffold_smiles="c1ccc2c(c1)CCN2",
                fragments=["C", "CC", "CCC"],
                patent_ids=["PAT003"],
            ),
        }

        # Sample core records
        self._cores = [
            CoreRecord(
                core_id="CORE001",
                core_smiles="c1ccccc1",
                connection_sites=6,
                patent_ids=["PAT001", "PAT002"],
            ),
            CoreRecord(
                core_id="CORE002",
                core_smiles="c1ccc2c(c1)cccc2",  # Naphthalene
                connection_sites=8,
                patent_ids=["PAT004"],
            ),
            CoreRecord(
                core_id="CORE003",
                core_smiles="c1ccncc1",  # Pyridine
                connection_sites=5,
                patent_ids=["PAT005"],
            ),
            CoreRecord(
                core_id="CORE004",
                core_smiles="c1ccc2[nH]ccc2c1",  # Indole
                connection_sites=7,
                patent_ids=["PAT006"],
            ),
            CoreRecord(
                core_id="CORE005",
                core_smiles="c1ccc2oc(=O)ccc2c1",  # Coumarin
                connection_sites=6,
                patent_ids=["PAT007"],
            ),
        ]

        # Fragment library for fallback generation
        self._fragment_library = [
            "C",
            "CC",
            "CCC",
            "F",
            "Cl",
            "Br",
            "N",
            "NC",
            "NCC",
            "O",
            "OC",
            "OCC",
            "NC(=O)C",
            "NC(=O)CC",
            "c1ccccc1",
            "c1ccncc1",
            "C(=O)O",
            "C(=O)N",
            "S",
            "SC",
            "CF",
            "C(F)(F)F",
            "c1ccc(O)cc1",
            "c1ccc(N)cc1",
        ]

    def search_by_fragment(self, fragment_smiles: str) -> list[PatentRecord]:
        """Search for patents containing the given fragment."""
        results = []
        for patent in self._patents.values():
            if patent.fragments and fragment_smiles in patent.fragments:
                results.append(patent)
        return results

    def search_by_scaffold(self, scaffold_smiles: str) -> list[ScaffoldRecord]:
        """Search for scaffolds matching the given SMILES."""
        if scaffold_smiles in self._scaffolds:
            return [self._scaffolds[scaffold_smiles]]
        return []

    def get_alternative_fragments(self, scaffold_smiles: str) -> list[str]:
        """Get alternative fragments for a scaffold."""
        if scaffold_smiles in self._scaffolds:
            return self._scaffolds[scaffold_smiles].fragments
        return []

    def search_cores(self, min_sites: int = 1) -> list[CoreRecord]:
        """Search for cores with minimum connection sites."""
        return [core for core in self._cores if core.connection_sites >= min_sites]

    def get_alternative_cores(self, core_smiles: str) -> list[CoreRecord]:
        """Get alternative cores."""
        # For mock, return all cores except the one provided
        return [core for core in self._cores if core.core_smiles != core_smiles]

    def get_fragment_library(self) -> list[str]:
        """Get the full fragment library for candidate generation."""
        return self._fragment_library.copy()

    def add_fragment(self, smiles: str) -> None:
        """Add a fragment to the library."""
        if smiles not in self._fragment_library:
            self._fragment_library.append(smiles)

    def add_core(self, core_smiles: str, connection_sites: int, patent_ids: Optional[list[str]] = None) -> None:
        """Add a core to the database."""
        core_id = f"CORE{len(self._cores) + 1:03d}"
        self._cores.append(
            CoreRecord(
                core_id=core_id,
                core_smiles=core_smiles,
                connection_sites=connection_sites,
                patent_ids=patent_ids or [],
            )
        )
