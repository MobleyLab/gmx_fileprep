from openff.toolkit import Molecule


def create_molecule(mol, ligcode="LIG"):
    ligand_molecule = Molecule.from_rdkit(mol)
    ligand_molecule.name = ligcode
    for atom in ligand_molecule.atoms:
        atom.metadata["residue_name"] = ligcode
    return ligand_molecule


def create_molecule(mol, ligcode="LIG"):
    ligand_molecule = Molecule.from_mol2(mol)
    ligand_molecule.name = ligcode
    for atom in ligand_molecule.atoms:
        atom.metadata["residue_name"] = ligcode
    return ligand_molecule