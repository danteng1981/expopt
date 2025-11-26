"""RGroup splicing implementation for attaching fragments to scaffolds."""

from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions


class RGroupSplicing:
    """
    Attach replacement fragments to scaffolds via single bonds.

    Handles attachment at specified positions marked with dummy atoms
    or at specified atom indices.
    """

    def __init__(self):
        """Initialize RGroup splicing utility."""
        pass

    def splice_fragment(
        self,
        scaffold_smiles: str,
        fragment_smiles: str,
        scaffold_attach_idx: Optional[int] = None,
        fragment_attach_idx: Optional[int] = None,
    ) -> Optional[str]:
        """
        Attach a fragment to a scaffold via single bond.

        Args:
            scaffold_smiles: SMILES of scaffold (may contain attachment points [*])
            fragment_smiles: SMILES of fragment to attach
            scaffold_attach_idx: Atom index on scaffold for attachment (if no [*])
            fragment_attach_idx: Atom index on fragment for attachment (if no [*])

        Returns:
            SMILES of combined molecule or None if splicing fails
        """
        scaffold_mol = Chem.MolFromSmiles(scaffold_smiles)
        fragment_mol = Chem.MolFromSmiles(fragment_smiles)

        if scaffold_mol is None or fragment_mol is None:
            return None

        # Look for dummy atoms in scaffold
        scaffold_attach = self._find_attachment_point(scaffold_mol)
        if scaffold_attach is None:
            scaffold_attach = scaffold_attach_idx

        # Look for dummy atoms in fragment
        fragment_attach = self._find_attachment_point(fragment_mol)
        if fragment_attach is None:
            fragment_attach = fragment_attach_idx

        # Use first available position if still None
        if scaffold_attach is None:
            scaffold_attach = 0
        if fragment_attach is None:
            fragment_attach = 0

        return self._connect_molecules(
            scaffold_mol, fragment_mol, scaffold_attach, fragment_attach
        )

    def _find_attachment_point(self, mol: Chem.Mol) -> Optional[int]:
        """
        Find attachment point (dummy atom [*]) in molecule.

        Args:
            mol: RDKit molecule

        Returns:
            Atom index of attachment point or None
        """
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:  # Dummy atom
                # Return the neighboring atom index
                neighbors = atom.GetNeighbors()
                if neighbors:
                    return neighbors[0].GetIdx()
        return None

    def _connect_molecules(
        self,
        mol1: Chem.Mol,
        mol2: Chem.Mol,
        attach1: int,
        attach2: int,
    ) -> Optional[str]:
        """
        Connect two molecules via single bond.

        Args:
            mol1: First molecule
            mol2: Second molecule
            attach1: Attachment point on first molecule
            attach2: Attachment point on second molecule

        Returns:
            SMILES of combined molecule
        """
        # Remove dummy atoms from both molecules
        mol1_clean = self._remove_dummy_atoms(mol1)
        mol2_clean = self._remove_dummy_atoms(mol2)

        if mol1_clean is None or mol2_clean is None:
            return None

        # Combine molecules
        combined = Chem.CombineMols(mol1_clean, mol2_clean)
        combined_edit = Chem.RWMol(combined)

        # Add bond between attachment points
        # Adjust attach2 index for the combined molecule
        attach2_adjusted = attach2 + mol1_clean.GetNumAtoms()

        try:
            combined_edit.AddBond(attach1, attach2_adjusted, Chem.BondType.SINGLE)
            result = combined_edit.GetMol()
            Chem.SanitizeMol(result)
            return Chem.MolToSmiles(result)
        except Exception:
            return None

    def _remove_dummy_atoms(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """
        Remove dummy atoms from molecule.

        Args:
            mol: RDKit molecule

        Returns:
            Molecule without dummy atoms
        """
        mol_copy = Chem.RWMol(mol)

        # Find dummy atoms
        dummy_indices = []
        for atom in mol_copy.GetAtoms():
            if atom.GetAtomicNum() == 0:
                dummy_indices.append(atom.GetIdx())

        # Remove in reverse order to maintain indices
        for idx in sorted(dummy_indices, reverse=True):
            mol_copy.RemoveAtom(idx)

        try:
            result = mol_copy.GetMol()
            Chem.SanitizeMol(result)
            return result
        except Exception:
            return mol  # Return original if sanitization fails

    def splice_with_reaction(
        self,
        scaffold_smiles: str,
        fragment_smiles: str,
    ) -> list[str]:
        """
        Splice using reaction SMARTS for more control.

        Creates a generic attachment reaction that joins
        molecules at marked positions.

        Args:
            scaffold_smiles: Scaffold with attachment point
            fragment_smiles: Fragment with attachment point

        Returns:
            List of product SMILES
        """
        # Generic single bond formation reaction
        # [*:1]-[!$([*])] reacts with [*:2]-[!$([*])]
        rxn_smarts = "[*:1][*:2].[*:3][*:4]>>[*:1][*:4]"

        try:
            rxn = rdChemReactions.ReactionFromSmarts(rxn_smarts)
            scaffold_mol = Chem.MolFromSmiles(scaffold_smiles)
            fragment_mol = Chem.MolFromSmiles(fragment_smiles)

            if scaffold_mol is None or fragment_mol is None:
                return []

            products = rxn.RunReactants((scaffold_mol, fragment_mol))

            result_smiles = []
            for product_set in products:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        smiles = Chem.MolToSmiles(product)
                        if smiles not in result_smiles:
                            result_smiles.append(smiles)
                    except Exception:
                        continue

            return result_smiles
        except Exception:
            return []

    def enumerate_attachments(
        self,
        scaffold_smiles: str,
        fragment_smiles: str,
    ) -> list[str]:
        """
        Enumerate all possible attachment configurations.

        Args:
            scaffold_smiles: Scaffold SMILES
            fragment_smiles: Fragment SMILES

        Returns:
            List of unique product SMILES
        """
        scaffold_mol = Chem.MolFromSmiles(scaffold_smiles)
        fragment_mol = Chem.MolFromSmiles(fragment_smiles)

        if scaffold_mol is None or fragment_mol is None:
            return []

        products = []
        seen_smiles = set()

        # Try all combinations of attachment points
        for s_idx in range(scaffold_mol.GetNumAtoms()):
            for f_idx in range(fragment_mol.GetNumAtoms()):
                result = self.splice_fragment(
                    scaffold_smiles, fragment_smiles, s_idx, f_idx
                )
                if result and result not in seen_smiles:
                    products.append(result)
                    seen_smiles.add(result)

        return products
