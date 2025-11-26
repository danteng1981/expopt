"""Constraint checking utilities for fragment/core selection."""

from rdkit import Chem


class AtomCountConstraint:
    """Check atom count constraints for fragments and cores."""

    def __init__(self, target_atom_count: int, min_offset: int = 3, max_multiplier: float = 1.5):
        """
        Initialize constraint checker.

        Args:
            target_atom_count: Number of atoms in target fragment/core
            min_offset: Minimum offset below target (default: 3)
            max_multiplier: Maximum multiplier for upper bound (default: 1.5)
        """
        self.target_atom_count = target_atom_count
        self.min_count = max(1, target_atom_count - min_offset)
        self.max_count = int(target_atom_count * max_multiplier)

    def check_constraint(self, mol: Chem.Mol) -> bool:
        """
        Check if molecule meets atom count constraint.

        Args:
            mol: RDKit molecule object

        Returns:
            True if constraint is met, False otherwise
        """
        if mol is None:
            return False

        atom_count = mol.GetNumHeavyAtoms()
        return self.min_count <= atom_count <= self.max_count

    def filter_candidates(self, candidates: list[Chem.Mol]) -> list[Chem.Mol]:
        """
        Filter candidates by atom count constraint.

        Args:
            candidates: List of candidate molecules

        Returns:
            List of molecules meeting the constraint
        """
        return [mol for mol in candidates if self.check_constraint(mol)]


class ConnectionSiteConstraint:
    """Check connection site constraints for core hopping."""

    def __init__(self, min_sites: int):
        """
        Initialize connection site constraint.

        Args:
            min_sites: Minimum number of potential connection sites required
        """
        self.min_sites = min_sites

    def count_connection_sites(self, mol: Chem.Mol) -> int:
        """
        Count potential connection sites in a molecule.

        Connection sites are identified as atoms with attachment points
        (dummy atoms [*] or atoms with explicit single bond attachments).

        Args:
            mol: RDKit molecule object

        Returns:
            Number of potential connection sites
        """
        if mol is None:
            return 0

        sites = 0

        for atom in mol.GetAtoms():
            # Check for dummy atoms (attachment points)
            if atom.GetAtomicNum() == 0:
                sites += 1
            # Check for atoms that could accept additional bonds
            elif self._can_accept_bond(atom):
                sites += 1

        return sites

    def _can_accept_bond(self, atom: Chem.Atom) -> bool:
        """
        Check if an atom can accept additional bonds.

        Args:
            atom: RDKit atom object

        Returns:
            True if atom can accept another bond
        """
        # Get current valence
        total_valence = atom.GetTotalValence()
        default_valence = Chem.GetPeriodicTable().GetDefaultValence(atom.GetAtomicNum())

        # Handle multiple possible valences
        if isinstance(default_valence, tuple):
            max_valence = max(default_valence)
        else:
            max_valence = default_valence

        return total_valence < max_valence

    def check_constraint(self, mol: Chem.Mol) -> bool:
        """
        Check if molecule meets connection site constraint.

        Args:
            mol: RDKit molecule object

        Returns:
            True if constraint is met, False otherwise
        """
        return self.count_connection_sites(mol) >= self.min_sites

    def filter_candidates(self, candidates: list[Chem.Mol]) -> list[Chem.Mol]:
        """
        Filter candidates by connection site constraint.

        Args:
            candidates: List of candidate molecules

        Returns:
            List of molecules meeting the constraint
        """
        return [mol for mol in candidates if self.check_constraint(mol)]
