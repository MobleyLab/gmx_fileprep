from openff.toolkit import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField

def setup_ff(offxml):
    return ForceField(offxml)


def create_molecule(mol2, ligcode="LIG"):
    ligand_molecule = Molecule(mol2)
    ligand_molecule.name = ligcode
    for atom in ligand_molecule.atoms:
        atom.metadata["residue_name"] = ligcode
    return ligand_molecule