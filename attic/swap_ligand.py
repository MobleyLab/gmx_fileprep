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
import tempfile
import time
from pathlib import Path
import parmed
import os
from rdkit import Chem
from rdkit.Chem import AllChem

def main():
    start_time = time.time()

    ###################
    ## COMPLEX SETUP ##
    ###################

    # print("\nComplex Setup\n" + "="*13)
    # print("Step  1: Loading parameterized Parmed Structure for Complex + Solvent + Ions")
    # pmd_receptor_struct = parmed.load_file(complex_top, xyz=complex_gro)

    # print("Step  2: Creating a ligand .mol2 file using Complex LIG positions")
    # complex_lig_struct = pmd_receptor_struct[pmdw.create_mask_str([complex_ligcode])]
    # with tempfile.NamedTemporaryFile(suffix=".mol2") as complex_mol2:
    #     complex_lig_struct.save(complex_mol2.name, overwrite=True)
    #     pmdw.edit_mol2_positions(ligand_mol2, complex_mol2.name, ligand_mol2_complex)

    # print("Step  3: Parameterizing Ligand with OpenFF")
    # openff_ff = offw.setup_ff(ligand_ff)
    # lig_mol = offw.create_molecule(ligand_mol2_complex, complex_ligcode)
    # lig_positions = lig_mol.conformers[0]
    # lig_topology = lig_mol.to_topology()
    # lig_system = openff_ff.create_openmm_system(lig_topology)

    # print("Step  4: Creating Parmed Structure for Ligand")
    # pmd_lig_struct = pmdw.create_pmd_ligand(lig_topology, lig_system, lig_positions)

    # print("Step  5: Creating Parmed Structure for Solvated Complex: Protein + Ligand + Solvent + Ions")
    # exclude_resids = [complex_ligcode, 'Na+', 'Cl-', 'HOH']
    # unq_lig_struct = pmdw.unique_atom_types(pmd_lig_struct, lig_mol.name)
    # unq_lig_struct = pmdw.fix_gen_pairs(unq_lig_struct)
    # unq_lig_struct.save("/Users/megosato/Desktop/test.top", overwrite=True)

    # pmd_receptor_struct.save("/Users/megosato/Desktop/test2.top", overwrite=True)

    # mask = pmdw.create_mask_from_exclusion(pmd_receptor_struct, exclude_resids)
    # pmd_complex_struct = pmd_receptor_struct[mask] \
    #                         + unq_lig_struct \
    #                         + pmd_receptor_struct[':Na+:Cl-:HOH']

    # # For some reason... the parmed complex needs to be saved out before running
    # # unique_atom_types() in order to avoid the scaling 1-4 error
    # tmp_top = tempfile.NamedTemporaryFile(suffix='.top', delete=False)
    # tmp_gro = tempfile.NamedTemporaryFile(suffix='.gro', delete=False)
    # pmd_complex_struct.save(tmp_top.name, overwrite=True)
    # pmd_complex_struct.save(tmp_gro.name, overwrite=True)
    # pmd_complex_struct = parmed.load_file(tmp_top.name, xyz=tmp_gro.name)
    # os.unlink(tmp_gro.name)
    # os.unlink(tmp_top.name)

    # print("Step  6: Ensuring All Atomtypes are Unique")
    # unq_complex_struct = pmdw.unique_atom_types(pmd_complex_struct, lig_mol.name)
    # final_complex_struct = pmdw.fix_gen_pairs(unq_complex_struct)

    # print("Step  7: Saving out Complex .gro and .top and .pdb")
    # final_complex_struct.save(complex_gro_out, overwrite=True)
    # final_complex_struct.save(complex_top_out, overwrite=True)
    # final_complex_struct.save(complex_pdb_out, overwrite=True)

    # # Add a sanity check that confirms gro files are the same before and after addition
    # # of differently parameterized ligand
    # print("Step  8: Sanity check:")
    # print("  .gro: modified file should columns be a subset of the original (subset or !subset)")
    # print("  .top: modified file should be different from original")
    # print(f"    filetype  original  modified  difference")
    # print(f"    ========  ========  ========  ==========")
    # print(f"      .gro    {iow.count_lines(complex_gro):>8}  {iow.count_lines(complex_gro_out):>8}  {'subset' if iow.is_a_subset(complex_gro_out, complex_gro) else '!subset':>10}")
    # print(f"      .top    {iow.count_lines(complex_top):>8}  {iow.count_lines(complex_top_out):>8}  {iow.diff(complex_top, complex_top_out):>10}")



    ###################
    ## SOLVENT SETUP ##
    ###################
    print("\n\nSolvent Setup\n" + "="*13)

    print("Step  1: Loading parameterized Parmed Structure for Ligand + Solvent + Ions")
    pmd_solvent_struct = parmed.load_file(solvent_top, xyz=solvent_gro)



    print("Step  2: Creating a ligand .mol2 file using Solvent UNL positions")

    rdmol = None
    solvent_lig_struct = pmd_solvent_struct[pmdw.create_mask_str([solvent_ligcode])]
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
    openff_ff = offw.setup_ff(ligand_ff)
    lig_mol = offw.create_molecule(rdmol, solvent_ligcode)
    lig_positions = lig_mol.conformers[0]
    lig_topology = lig_mol.to_topology()
    lig_system = openff_ff.create_openmm_system(lig_topology)

    print("Step  4: Creating Parmed Structure for Ligand")
    pmd_lig_struct = pmdw.create_pmd_ligand(lig_topology, lig_system, lig_positions)
    #pmd_lig_struct.save("/Users/megosato/Desktop/lig_after_ff.top")

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
    final_solvent_struct.save(solvent_gro_out, overwrite=True)
    final_solvent_struct.save(solvent_top_out, overwrite=True)


    # Add a sanity check that confirms gro files are the same before and after addition
    # of differently parameterized ligand
    print("Step  9: Sanity check:")
    print("  .gro: modified file columns should be a subset of the original (subset or !subset)")
    print("  .top: modified file should be different from original")
    print(f"    filetype  original  modified  difference")
    print(f"    ========  ========  ========  ==========")
    print(f"      .gro    {iow.count_lines(solvent_gro):>8}  {iow.count_lines(solvent_gro_out):>8}  {'subset' if iow.is_a_subset(solvent_gro_out, solvent_gro) else '!subset':>10}")
    print(f"      .top    {iow.count_lines(solvent_top):>8}  {iow.count_lines(solvent_top_out):>8}  {iow.diff(solvent_top, solvent_top_out):>10}")


    # end_time = time.time()
    # print(f"Total Time Elapsed: {(end_time - start_time)/60:.2} min")

if __name__ == "__main__":

    desktop_path = Path("/Users/megosato/Desktop")
    si_path = desktop_path / "SI/BACE1/input"
    output_path = desktop_path / "lig11_bace_prep"

    ligand_ff = desktop_path / "bace_prep/sage-2.1.0rc.offxml"

    complex_ligcode = "LIG"
    solvent_ligcode = "UNL"

    for ligfam in ["amide_series", "spirocycles", "biphenyl"]:
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

            if lig not in ["lig_41", "lig_67"]:
                continue

            lig_si_dir = ligfam_dir / lig
            complex_gro = str(lig_si_dir / "complex.gro")
            complex_top = str(lig_si_dir / "complex.top")
            solvent_gro = str(lig_si_dir / "solvent.gro")
            solvent_top = str(lig_si_dir / "solvent.top")
            ligand_mol2 = str(mol2_dir / f"{lig}.mol2")

            if ligfam == "biphenyl":
                lig = "bJ_" + lig.split("_")[-1]

            lig_output_dir = fam_output_dir / lig
            lig_output_dir.mkdir(parents=True, exist_ok=True)

            ligand_mol2_complex = str(lig_output_dir / f"{lig}_complex.mol2")
            ligand_mol2_solvent = str(lig_output_dir / f"{lig}_solvent.mol2")
            complex_gro_out = str(lig_output_dir / "complex.gro")
            complex_top_out = str(lig_output_dir / "complex.top")
            complex_pdb_out = str(lig_output_dir / "complex.pdb")
            solvent_gro_out = str(lig_output_dir / "solvent.gro")
            solvent_top_out = str(lig_output_dir / "solvent.top")

            main()