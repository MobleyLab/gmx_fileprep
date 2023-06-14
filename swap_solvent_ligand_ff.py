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

def main(ligand_ff, gro, top, ligcode, ligand_mol2, gro_out, top_out):
    start_time = time.time()

    ###################
    ## SOLVENT SETUP ##
    ###################
    print("\n\nSolvent Setup\n" + "="*13)

    print("Step  1: Loading parameterized Parmed Structure for Ligand + Solvent + Ions")
    pmd_solvent_struct = parmed.load_file(top, xyz=gro)

    print("Step  2: Creating a ligand .mol2 file using Solvent UNL positions")
    utils.printerr("    WARNING: NEED TO FIX the Hs under a different resid than the rest of the ligand")
    utils.printerr("    WARNING: NEED TO FIX atomtypes have an x (i.e. C1x)")

    rdmol = None
    solvent_lig_struct = pmd_solvent_struct[pmdw.create_mask_str([ligcode])]
    with tempfile.NamedTemporaryFile(suffix=".pdb") as solvent_pdb:
        solvent_lig_struct.save(solvent_pdb.name, overwrite=True)
        #pmdw.edit_mol2_positions(ligand_mol2, solvent_mol2.name, ligand_mol2_solvent)

        mol2 = Chem.MolFromMol2File(ligand_mol2)
        mol2_smiles = Chem.MolToSmiles(mol2)
        template = Chem.MolFromSmiles(mol2_smiles)
        docked_pose = AllChem.MolFromPDBFile(solvent_pdb.name)

        #Assign the bond order to force correct valence
        rdmol = AllChem.AssignBondOrdersFromTemplate(template, docked_pose)

    print("Step  3: Parameterizing Ligand with OpenFF")
    openff_ff = ForceField(ligand_ff)
    lig_mol = offw.create_molecule(rdmol, ligcode)
    lig_positions = lig_mol.conformers[0]
    lig_topology = lig_mol.to_topology()
    lig_system = openff_ff.create_openmm_system(lig_topology)

    print("Step  4: Creating Parmed Structure for Ligand")
    pmd_lig_struct = pmdw.create_pmd_ligand(lig_topology, lig_system, lig_positions)
    pmd_lig_struct = pmdw.remove_x_atomname(pmd_lig_struct)

    print("Step  5: Creating Parmed Structure for Solvent: Ligand + Solvent + Ions")

    pmd_solvent_struct = pmd_solvent_struct[':Cl-:Na+'] \
                        + pmd_lig_struct \
                        + pmd_solvent_struct[':HOH']

    unq_solvent_struct = pmdw.unique_atom_types(pmd_solvent_struct, lig_mol.name)
    #unq_solvent_struct.save("/Users/megosato/Desktop/solvent_after_combo.top")

    print("Step  6: Ensuring Ligand Atomtypes are Unique")
    #unq_solvent_struct = pmdw.unique_atom_types(pmd_solvent_struct, lig_mol.name)
    final_solvent_struct = pmdw.fix_gen_pairs(unq_solvent_struct)

    print("Step  7: Saving out Solvent .gro and .top")
    final_solvent_struct.save(gro_out, overwrite=True)
    final_solvent_struct.save(top_out, overwrite=True)


    # Add a sanity check that confirms gro files are the same before and after addition
    # of differently parameterized ligand
    print("Step  9: Sanity check:")
    print("  .gro: modified file columns should be a subset of the original (subset or !subset)")
    print("  .top: modified file should be different from original")
    print(f"    filetype  original  modified  difference")
    print(f"    ========  ========  ========  ==========")
    print(f"      .gro    {iow.count_lines(gro):>8}  {iow.count_lines(gro_out):>8}  {'subset' if iow.is_a_subset(gro_out, gro) else '!subset':>10}")
    print(f"      .top    {iow.count_lines(top):>8}  {iow.count_lines(top_out):>8}  {iow.diff(top, top_out):>10}")


    # end_time = time.time()
    # print(f"Total Time Elapsed: {(end_time - start_time)/60:.2} min")

if __name__ == "__main__":

    ligand_ff = "openff_unconstrained-2.0.0.offxml"
    ligcode = "LIG"

    desktop_path = Path("/Users/megosato/Desktop")
    si_path = desktop_path / "SI/BACE1/input"

    ligfam = "biphenyl"
    ligfam_dir = si_path / ligfam
    mol2_dir = ligfam_dir / "mol2"

    lig_names = ["lig_04", "lig_02", "lig_03"]

    fam_output_dir = desktop_path / "biphenyl_2.0.0_solv_redo"
    fam_output_dir.mkdir(parents=True, exist_ok=True)

    for lig in lig_names:
        print(ligfam, lig)
        lig_si_dir = ligfam_dir / lig
        gro = str(lig_si_dir / "solvent.gro")
        top = str(lig_si_dir / "solvent.top")
        ligand_mol2 = str(mol2_dir / f"{lig}.mol2")

        lig_output_dir = fam_output_dir / lig
        lig_output_dir.mkdir(parents=True, exist_ok=True)

        gro_out = str(lig_output_dir / "solvent.gro")
        top_out = str(lig_output_dir / "solvent.top")
        pdb_out = str(lig_output_dir / "solvent.pdb")

        main(ligand_ff, gro, top, ligcode, ligand_mol2, gro_out, top_out)

    # desktop_path = Path("/Users/megosato/Desktop")
    # si_path = desktop_path / "SI/BACE1/input"
    # output_path = desktop_path / "TESTING"

    # ligand_ff = desktop_path / "bace_prep/sage-2.1.0rc.offxml"

    # ligcode = "UNL"

    # for ligfam in ["amide_series"]:#, "spirocycles", "biphenyl"]:
    #     ligfam_dir = si_path / ligfam
    #     mol2_dir = ligfam_dir / "mol2"
    #     lig_names = []
    #     mol2_files = Path(mol2_dir).glob('*')
    #     for f in mol2_files:
    #         lig_names.append(str(f.stem))

    #     fam_output_dir = output_path / ligfam
    #     fam_output_dir.mkdir(parents=True, exist_ok=True)

    #     for lig in lig_names:

    #         print(ligfam, lig)

    #         lig_si_dir = ligfam_dir / lig
    #         gro = str(lig_si_dir / "solvent.gro")
    #         top = str(lig_si_dir / "solvent.top")
    #         ligand_mol2 = str(mol2_dir / f"{lig}.mol2")

    #         lig_output_dir = fam_output_dir / lig
    #         lig_output_dir.mkdir(parents=True, exist_ok=True)

    #         gro_out = str(lig_output_dir / "solvent.gro")
    #         top_out = str(lig_output_dir / "solvent.top")

    #         main(ligand_ff, gro, top, ligcode, ligand_mol2, gro_out, top_out)

