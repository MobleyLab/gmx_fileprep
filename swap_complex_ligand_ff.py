# swap_complex_ligand_ff.py
#  * this is an EXAMPLE main file using the wrapper functions
#  * workflow creates ligand and complex .gro and .top files
#  * complex files are ready for energy minimization
#  * ligand files need to be solvated and have ions added
#    the water and ion FF parameters have already been added
#  * the default parameters in the gromacs .top file are set 
#    for SepTop (gen-pairs: no, fudgeLJ: 0.5)
#
# THIS IS A WIP... some steps still are not fully functional

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

def main(ligand_ff, gro, top, ligcode, ligand_mol2, gro_out, top_out, pdb_out):
    '''
        ligand_ff:   ligand forcefield file path or filename
            * use ligand ff file path if the ff has not been released
            * use filename if ff is available through openff package
        gro:         .gro file of complex with old ff
        top:         .top file of complex with old ff
        ligcode:     3 letter code for ligand
        ligand_mol2: .mol2 file for ligand
        gro_out:     .gro file to write newly parameterized complex 
        top_out:     .top file to write newly parameterized complex
        pdb_out:     .pdb file to write newly parameterized complex
    '''

    start_time = time.time()

    ##################
    # COMPLEX SETUP ##
    ##################

    print("\nComplex Setup\n" + "="*13)
    print("Step  1: Loading parameterized Parmed Structure for Complex + Solvent + Ions")
    pmd_receptor_struct = parmed.load_file(top, xyz=gro)

    print("Step  2: Creating a ligand .mol2 file using Complex LIG positions")
    utils.printerr("    WARNING: Check the valence and bond order of your ligands")
    complex_lig_struct = pmd_receptor_struct[pmdw.create_mask_str([ligcode])]
    rdmol = pmdw.reorder_mol_atoms(complex_lig_struct, ligand_mol2)
    utils.printerr("    WARNING: rdkit mol H coords may be slightly shifted from template")

    print("Step  3: Parameterizing Ligand with OpenFF")
    openff_ff = ForceField(ligand_ff)
    lig_mol = Molecule.from_rdkit(rdmol)
    lig_positions = lig_mol.conformers[0]
    lig_topology = lig_mol.to_topology()
    lig_system = openff_ff.create_openmm_system(lig_topology)

    print("Step  4: Creating Parmed Structure for Ligand")
    pmd_lig_struct = pmdw.create_pmd_ligand(lig_topology, lig_system, lig_positions)
    # do some sketchy stuff to fix atom names with x and Hs in a different residue
    pmd_lig_struct = pmdw.remove_x_atomname(pmd_lig_struct)

    with tempfile.NamedTemporaryFile(suffix='.mol2') as complex_mol2:
        complex_lig_struct.save(complex_mol2.name, overwrite=True)
        pmd_lig_struct = pmdw.edit_mol2_positions(pmd_lig_struct, complex_mol2.name) 

    print("Step  5: Creating Parmed Structure for Solvated Complex: Protein + Ligand + Solvent + Ions")
    exclude_resids = [ligcode, 'Na+', 'Cl-', 'HOH']
    unq_lig_struct = pmdw.unique_atom_types(pmd_lig_struct, lig_mol.name)
    unq_lig_struct = pmdw.fix_gen_pairs(unq_lig_struct)

    mask = pmdw.create_mask_from_exclusion(pmd_receptor_struct, exclude_resids)
    pmd_complex_struct = pmd_receptor_struct[mask] \
                            + unq_lig_struct \
                            + pmd_receptor_struct[':Na+:Cl-:HOH']

    # For some reason... the parmed complex needs to be saved out before running
    # unique_atom_types() in order to avoid the scaling 1-4 error
    tmp_top = tempfile.NamedTemporaryFile(suffix='.top', delete=False)
    tmp_gro = tempfile.NamedTemporaryFile(suffix='.gro', delete=False)
    pmd_complex_struct.save(tmp_top.name, overwrite=True)
    pmd_complex_struct.save(tmp_gro.name, overwrite=True)
    pmd_complex_struct = parmed.load_file(tmp_top.name, xyz=tmp_gro.name)
    os.unlink(tmp_gro.name)
    os.unlink(tmp_top.name)

    print("Step  6: Ensuring All Atomtypes are Unique")
    unq_complex_struct = pmdw.unique_atom_types(pmd_complex_struct, lig_mol.name)
    final_complex_struct = pmdw.fix_gen_pairs(unq_complex_struct)

    print("Step  7: Saving out Complex .gro and .top and .pdb")
    final_complex_struct.save(gro_out, overwrite=True)
    final_complex_struct.save(top_out, overwrite=True)
    final_complex_struct.save(pdb_out, overwrite=True)

    # Add a sanity check that confirms gro files are the same before and after addition
    # of differently parameterized ligand
    print("Step  8: EXTREMELY BASIC not comprehensive sanity check:")
    print("  .gro: modified file should columns be a subset of the original (subset or !subset)")
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
    gro = base_path / "SI/BACE1/input/biphenyl/lig_04/complex.gro"
    top = base_path / "/Users/megosato/Desktop/SI/BACE1/input/biphenyl/lig_04/complex.top"
    ligcode = "LIG"
    ligand_mol2 = base_path / "SI/BACE1/input/biphenyl/mol2/lig_04.mol2"
    gro_out = base_path / "complex.gro"
    top_out = base_path / "complex.top"
    pdb_out = base_path / "complex.pdb"

    main(
        ligand_ff, 
        str(gro), str(top), ligcode, str(ligand_mol2), 
        str(gro_out), str(top_out), str(pdb_out),
    )






