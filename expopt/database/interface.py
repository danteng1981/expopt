"""Abstract database interface for patent and scaffold lookups."""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional

from rdkit import Chem


@dataclass
class PatentRecord:
    """Represents a patent record containing molecular data."""

    patent_id: str
    smiles: str
    scaffold_smiles: Optional[str] = None
    fragments: Optional[list[str]] = None
    metadata: Optional[dict] = None


@dataclass
class ScaffoldRecord:
    """Represents a scaffold record with associated fragments."""

    scaffold_id: str
    scaffold_smiles: str
    fragments: list[str]  # SMILES of associated fragments
    patent_ids: list[str]  # Associated patents


@dataclass
class CoreRecord:
    """Represents a core skeleton record."""

    core_id: str
    core_smiles: str
    connection_sites: int
    patent_ids: list[str]
    metadata: Optional[dict] = None


class DatabaseInterface(ABC):
    """Abstract interface for database operations."""

    @abstractmethod
    def search_by_fragment(self, fragment_smiles: str) -> list[PatentRecord]:
        """
        Search for patent records containing the given fragment.

        Args:
            fragment_smiles: SMILES string of fragment to search

        Returns:
            List of matching patent records
        """
        pass

    @abstractmethod
    def search_by_scaffold(self, scaffold_smiles: str) -> list[ScaffoldRecord]:
        """
        Search for scaffold records matching the given scaffold.

        Args:
            scaffold_smiles: SMILES string of scaffold to search

        Returns:
            List of matching scaffold records
        """
        pass

    @abstractmethod
    def get_alternative_fragments(self, scaffold_smiles: str) -> list[str]:
        """
        Get alternative fragments that appear with the same scaffold.

        Args:
            scaffold_smiles: SMILES string of scaffold

        Returns:
            List of SMILES strings for alternative fragments
        """
        pass

    @abstractmethod
    def search_cores(self, min_sites: int = 1) -> list[CoreRecord]:
        """
        Search for core scaffolds in the database.

        Args:
            min_sites: Minimum number of connection sites

        Returns:
            List of core records
        """
        pass

    @abstractmethod
    def get_alternative_cores(self, core_smiles: str) -> list[CoreRecord]:
        """
        Get alternative cores used in similar patent contexts.

        Args:
            core_smiles: SMILES string of target core

        Returns:
            List of alternative core records
        """
        pass

    def smiles_to_mol(self, smiles: str) -> Optional[Chem.Mol]:
        """
        Convert SMILES to RDKit Mol object.

        Args:
            smiles: SMILES string

        Returns:
            RDKit Mol object or None if conversion fails
        """
        return Chem.MolFromSmiles(smiles)

    def mol_to_smiles(self, mol: Chem.Mol) -> Optional[str]:
        """
        Convert RDKit Mol to SMILES string.

        Args:
            mol: RDKit Mol object

        Returns:
            SMILES string or None if conversion fails
        """
        if mol is None:
            return None
        return Chem.MolToSmiles(mol)
