# swap_ligand_params.py
#  * this is an EXAMPLE main file using the wrapper functions
#  * workflow creates ligand and complex .gro and .top files
#  * complex files are ready for energy minimization
#  * ligand files need to be solvated and have ions added
#    the water and ion FF parameters have already been added
#  * the default parameters in the gromacs .top file are set 
#    for SepTop (gen-pairs: no, fudgeLJ: 0.5)

import parmed_wrapper as pmdw
import io_wrapper as iow
import openff_wrapper as offw
import utils

import tempfile
import time
from pathlib import Path
import os

import parmed
from rdkit import Chem
from rdkit.Chem import AllChem
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit import Molecule


def main(ligand_ff, gro, top, ligcode, ligand_mol2, gro_out, top_out):
    '''
        ligand_ff:   ligand forcefield file path or filename
            * use ligand ff file path if the ff has not been released
            * use filename if ff is available through openff package
        gro:         .gro file of solvent with old ff
        top:         .top file of solvent with old ff
        ligcode:     3 letter code for ligand
        ligand_mol2: .mol2 file for ligand
        gro_out:     .gro file to write newly parameterized solvent 
        top_out:     .top file to write newly parameterized solvent
    '''


    start_time = time.time()

    print("\n\nSolvent Setup\n" + "="*13)

    print("Step  1: Loading parameterized Parmed Structure for Ligand + Solvent + Ions")
    pmd_solvent_struct = parmed.load_file(top, xyz=gro)

    print("Step  2: Creating an rdkit mol of ligand using Solvent UNL positions and Mol2 atomtypes")
    solvent_lig_struct = pmd_solvent_struct[pmdw.create_mask_str([ligcode])]
    rdmol = pmdw.reorder_mol_atoms(solvent_lig_struct, ligand_mol2)
    utils.printerr("    WARNING: rdkit mol H coords may be slightly shifted from template")

    print("Step  3: Parameterizing Ligand with OpenFF")
    openff_ff = ForceField(ligand_ff)
    lig_mol = offw.create_molecule(rdmol, ligcode)
    lig_positions = lig_mol.conformers[0]
    lig_topology = lig_mol.to_topology()
    lig_system = openff_ff.create_openmm_system(lig_topology)
       
    print("Step  4: Creating Parmed Structure for Ligand")
    pmd_lig_struct = pmdw.create_pmd_ligand(lig_topology, lig_system, lig_positions)
    pmd_lig_struct = pmdw.remove_x_atomname(pmd_lig_struct)

    with tempfile.NamedTemporaryFile(suffix='.mol2') as solvent_mol2:
        solvent_lig_struct.save(solvent_mol2.name, overwrite=True)
        pmd_lig_struct = pmdw.edit_mol2_positions(pmd_lig_struct, solvent_mol2.name)


    print("Step  5: Creating Parmed Structure for Solvent: Ligand + Solvent + Ions")
    pmd_solvent_struct = pmd_solvent_struct[':Cl-:Na+'] \
                        + pmd_lig_struct \
                        + pmd_solvent_struct[':HOH']

    unq_solvent_struct = pmdw.unique_atom_types(pmd_solvent_struct, lig_mol.name)

    print("Step  6: Ensuring Ligand Atomtypes are Unique")
    final_solvent_struct = pmdw.fix_gen_pairs(unq_solvent_struct)

    print("Step  7: Saving out Solvent .gro and .top")
    final_solvent_struct.save(gro_out, overwrite=True)
    final_solvent_struct.save(top_out, overwrite=True)

    # Add a sanity check that confirms gro files are the same before and after addition
    # of differently parameterized ligand
    print("Step  9: Sanity check (VERY BASIC, PLEASE CHECK YOUR FILES):")
    print("  .gro: modified file columns should be a subset of the original (subset or !subset)")
    print("  .top: modified file should be different from original")
    print(f"    filetype  original  modified  difference")
    print(f"    ========  ========  ========  ==========")
    print(f"      .gro    {iow.count_lines(gro):>8}  {iow.count_lines(gro_out):>8}  {'subset' if iow.is_a_subset(gro_out, gro) else '!subset':>10}")
    print(f"      .top    {iow.count_lines(top):>8}  {iow.count_lines(top_out):>8}  {iow.diff(top, top_out):>10}")

    end_time = time.time()
    print(f"Total Time Elapsed: {(end_time - start_time)/60:.2} min")

if __name__ == "__main__":
    ligand_ff = "openff_unconstrained-2.1.0.offxml"

    base_path = Path("/Users/megosato/Desktop/")
    gro = base_path / "SI/BACE1/input/biphenyl/lig_04/solvent.gro"
    top = base_path / "/Users/megosato/Desktop/SI/BACE1/input/biphenyl/lig_04/solvent.top"
    ligcode = "UNL"
    ligand_mol2 = base_path / "SI/BACE1/input/biphenyl/mol2/lig_04.mol2"
    gro_out = base_path / "solvent.gro"
    top_out = base_path / "solvent.top"

    
    main(
        ligand_ff, 
        str(gro), str(top), 
        ligcode, 
        str(ligand_mol2), 
        str(gro_out), str(top_out),
    )

