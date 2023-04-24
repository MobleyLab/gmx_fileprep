# prep_complex.py
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

def main():
    start_time = time.time()
    print("Step  1: Creating parameterized Parmed Structure for Protein + Solvent + Ions")
    pmd_receptor_struct = pmdw.create_pmd_receptor(protein_pdb, protein_ff, water_ff)


    print("Step  2: Parameterizing Ligand with OpenFF")
    openff_ff = offw.setup_ff(ligand_ff)
    lig_mol = offw.create_molecule(ligand_mol2, ligcode)
    lig_positions = lig_mol.conformers[0]
    lig_topology = lig_mol.to_topology()
    lig_system = openff_ff.create_openmm_system(lig_topology)

    print("Step  3: Creating Parmed Structure for Ligand")
    pmd_lig_struct = pmdw.create_pmd_ligand(lig_topology, lig_system, lig_positions)

    print("Step  4: Ensuring Ligand Atomtypes are Unique")
    unq_lig_struct = pmdw.unique_atom_types(pmd_lig_struct, lig_mol.name)

    print("Step  5: Creating Parmed Structure for Solvated Complex: Protein + Ligand + Solvent + Ions")
    pmd_complex_struct = pmd_receptor_struct + pmd_lig_struct


    print("Step  6: Removing Clashing Waters")
    clashes = pmdw.find_clashing_water(pmd_complex_struct, lig_mol.name, 0.1)
    if len(clashes) != 0:
        clash_residues_str = ','.join([str(i) for i in clashes])
        print(f'    Removing ligand-clashing water residues {clash_residues_str}')
        pmd_complex_struct.strip(f':{clash_residues_str}')
    else:
        print("    No ligand-water clashes to resolve")

    print("Step  7: Renaming Residues and Atoms for compatibility with GMX commands")
    print("    This step may take a while (+5 minutes)...")
    unq_complex_struct = pmdw.unique_atom_types(pmd_complex_struct, lig_mol.name)
    sol_complex_struct = pmdw.rename_residue(unq_complex_struct, 'HOH', 'SOL')
    sol_atm_name_dict = {
        ('SOL', 'O'): 'OW',
        ('SOL', 'H1'): 'HW1',
        ('SOL', 'H2'): 'HW2',
    }
    sol_atm_complex_struct = pmdw.rename_atoms(sol_complex_struct, sol_atm_name_dict)
    final_complex_struct = pmdw.fix_gen_pairs(sol_atm_complex_struct)

    print("Step  8: Saving out Protein + Ligand + Solvent + Ions .gro and .top files")
    final_complex_struct.save(complex_gro_out, overwrite=True)
    final_complex_struct.save(complex_top_out, overwrite=True)

    print("Step  9: Creating parameterized Ligand toplogy from solvated complex topology")
    final_lig_struct = pmdw.complexsolv_to_ligsolv_top(
        final_complex_struct, 
        'LIG', 'SOL', 'NA', 'CL'
    )

    print("Step 10: Saving out Ligand .gro and .top files")
    unq_lig_struct.save(ligand_gro_out, overwrite=True)
    with tempfile.NamedTemporaryFile(suffix=".top") as tmp_top:
        final_lig_struct.save(tmp_top.name, overwrite=True)
        iow.rm_molecule(tmp_top.name, ligand_top_out, [ligcode])

    end_time = time.time()
    print(f"Total Time Elapsed: {(end_time - start_time)/60:.2} min")

if __name__ == "__main__":

    base = "/Users/megosato/Desktop/bace/"

    lig = "lig_16"

    ligcode = "UNL"

    protein_pdb = base + "6od6_spruce_chainA.pdb"
    protein_ff = 'amber14/protein.ff14SB.xml'
    water_ff = 'amber14/tip3p.xml'

    ligand_ff = base + "sage-2.1.0rc.offxml"
    ligand_mol2 = base + f"spirocycles/mol2/{lig}.mol2"

    complex_gro_out = base + f"spirocycles/{lig}/complex.gro"
    complex_top_out = base + f"spirocycles/{lig}/complex.top"

    ligand_gro_out = base + f"spirocycles/{lig}/ligand.gro"
    ligand_top_out = base + f"spirocycles/{lig}/ligand.top"

    main()






