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

def main(ligand_ff, gro, top, ligcode, ligand_mol2, gro_out, top_out, pdb_out):
    start_time = time.time()

    ##################
    # COMPLEX SETUP ##
    ##################

    print("\nComplex Setup\n" + "="*13)
    print("Step  1: Loading parameterized Parmed Structure for Complex + Solvent + Ions")
    pmd_receptor_struct = parmed.load_file(top, xyz=gro)

    print("Step  2: Creating a ligand .mol2 file using Complex LIG positions")
    utils.printerr("    WARNING: NEED TO FIX the Hs under a different resid than the rest of the ligand")
    utils.printerr("    WARNING: NEED TO FIX atomtypes have an x (i.e. C1x)")
    complex_lig_struct = pmd_receptor_struct[pmdw.create_mask_str([ligcode])]
    rdmol = None
    with tempfile.NamedTemporaryFile(suffix=".pdb") as complex_lig_pdb:
        complex_lig_struct.save(complex_lig_pdb.name, overwrite=True)

        mol2 = Chem.MolFromMol2File(ligand_mol2)
        mol2_smiles = Chem.MolToSmiles(mol2)
        template = Chem.MolFromSmiles(mol2_smiles)
        docked_pose = AllChem.MolFromPDBFile(complex_lig_pdb.name)

        #Assign the bond order to force correct valence
        rdmol = AllChem.AssignBondOrdersFromTemplate(template, docked_pose)

    print("Step  3: Parameterizing Ligand with OpenFF")
    openff_ff = ForceField(ligand_ff)
    lig_mol = offw.create_molecule(rdmol)
    lig_positions = lig_mol.conformers[0]
    lig_topology = lig_mol.to_topology()
    lig_system = openff_ff.create_openmm_system(lig_topology)

    print("Step  4: Creating Parmed Structure for Ligand")
    pmd_lig_struct = pmdw.create_pmd_ligand(lig_topology, lig_system, lig_positions)

    print("Step  5: Creating Parmed Structure for Solvated Complex: Protein + Ligand + Solvent + Ions")
    exclude_resids = [ligcode, 'Na+', 'Cl-', 'HOH']
    unq_lig_struct = pmdw.unique_atom_types(pmd_lig_struct, lig_mol.name)
    unq_lig_struct = pmdw.fix_gen_pairs(unq_lig_struct)
    unq_lig_struct.save("/Users/megosato/Desktop/test.top", overwrite=True)

    pmd_receptor_struct.save("/Users/megosato/Desktop/test2.top", overwrite=True)

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
    desktop_path = Path("/Users/megosato/Desktop")
    si_path = desktop_path / "SI/BACE1/input"
    output_path = desktop_path / "TESTING"

    ligand_ff = desktop_path / "bace_prep/sage-2.1.0rc.offxml"

    ligcode = "LIG"

    for ligfam in ["amide_series"]: #, "spirocycles", "biphenyl"]:
        ligfam_dir = si_path / ligfam
        mol2_dir = ligfam_dir / "mol2"
        lig_names = []
        mol2_files = Path(mol2_dir).glob('*')
        for f in mol2_files:
            lig_names.append(str(f.stem))

        fam_output_dir = output_path / ligfam
        fam_output_dir.mkdir(parents=True, exist_ok=True)

        for lig in lig_names:

            print(ligfam, lig)

            lig_si_dir = ligfam_dir / lig
            gro = str(lig_si_dir / "complex.gro")
            top = str(lig_si_dir / "complex.top")
            ligand_mol2 = str(mol2_dir / f"{lig}.mol2")

            lig_output_dir = fam_output_dir / lig
            lig_output_dir.mkdir(parents=True, exist_ok=True)

            gro_out = str(lig_output_dir / "complex.gro")
            top_out = str(lig_output_dir / "complex.top")
            pdb_out = str(lig_output_dir / "complex.pdb")

            main(ligand_ff, gro, top, ligcode, ligand_mol2, gro_out, top_out, pdb_out)