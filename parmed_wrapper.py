import os
import tempfile

import parmed
import mdtraj
from openmm import app
from openmm.unit import *

import utils

def create_pmd_receptor(protein_pdb, protein_ff, water_ff):
    pdbfile = app.PDBFile(str(protein_pdb))
    modeller = app.Modeller(pdbfile.topology, pdbfile.positions)

    ff = app.ForceField(protein_ff, water_ff)

    modeller.addSolvent(
        ff, 
        padding=1*nanometers, 
        model='tip3p', 
        ionicStrength=0.15*molar, 
        positiveIon='Na+', 
        negativeIon='Cl-'
    )

    protein_system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        rigidWater=False,
    )

    pmd_receptor_struct = parmed.openmm.load_topology(
        modeller.topology,
        protein_system,
        modeller.positions
    )
    return pmd_receptor_struct

def create_pmd_ligand(lig_topology, lig_system, lig_positions):
    return parmed.openmm.load_topology(
        lig_topology.to_openmm(),
        lig_system,
        lig_positions,
    )


def pmd_resid_key(resid):
    ''' converts a residue id string into a key that can be used
        to index a parmed structure (ex: 'HOH' --> ':HOH')
        
        Parameters:
        ===========
        resid:  str    residue id
        
        Returns:
        ========
        str     residue id as a key (prepended with ':')
    '''
    return f':{resid}'
    

def rename_residue(pmd_structure, resid, new_resid):
    ''' converts all molecules with the residue id, resid, in the 
        parmed structure to new name, new_resid
        
        Parameters:
        ===========
        pmd_structure:  parmed topology
        resid:          str              residue id to be changed
        new_resid:      str              new residue id to use
        
        Returns:
        ========
        parmed topology   new parmed topology object with the residue
                          ids renamed
    '''
    new_structure = pmd_structure
    for res in new_structure.residues:
        if res.name == resid:
            res.name = new_resid
    return new_structure

def rename_atoms(pmd_structure, rename_dict):
    ''' renames atoms in the parmed structure based on the rename_dict
        
        Parameters:
        ===========
        pmd_structure:  parmed topology
        rename_dict:    dict             {('resid', 'old atm name'): 'new atm name'}
        
        Returns:
        ========
        parmed topology   new parmed topology object with the atoms
                          renamed
    '''
    new_structure = pmd_structure
    for r in new_structure.residues:
        for a in r.atoms:
            for resid,atomfrom in rename_dict:
                if r.name == resid and a.name == atomfrom:
                    a.name = rename_dict[(resid,atomfrom)]
    return new_structure


def complexsolv_to_ligsolv_top(pmd_structure, ligresid, solresid, ion1resid, ion2resid):
    ''' create a structure of the parameterized ligand and solvent 
        molecules from the solvated complex topology without any of 
        the protein atomtypes or FF parameters
        
        Parameters:
        ===========
        pmd_structure:  parmed topology   solvated complex parmed topology object 
        ligresid:       str               3 letter ligand residue id
        solresid:       str               solvent residue id
        ion1resid:      str               ion1 residue id
        ion2resid:      str               ion2 residue id
        
        Returns:
        ========
        parmed topology   ligand + solvent parmed topology object with ligand, 
                          solvent, and ions FF params
        
    '''
    
    ligsol_structure = pmd_structure[pmd_resid_key(ligresid)] \
                        + pmd_structure[pmd_resid_key(solresid)] \
                        + pmd_structure[pmd_resid_key(ion1resid)] \
                        + pmd_structure[pmd_resid_key(ion2resid)]
    
    return ligsol_structure

def fix_gen_pairs(structure):
    ''' modifies the .top file so that gen-pairs and fudgeLJ options
        are yes and 0.5 respecively (required for SepTop)
    '''
    tmp_top = tempfile.NamedTemporaryFile(suffix='.top', delete=False)
    tmp_modtop = tempfile.NamedTemporaryFile(suffix='.top', delete=False)
    tmp_gro = tempfile.NamedTemporaryFile(suffix='.gro', delete=False)
    
    structure.save(tmp_top.name, overwrite=True)
    structure.save(tmp_gro.name, overwrite=True)

    tmp_top_f = open(tmp_top.name, 'r')
    tmp_modtop_f = open(tmp_modtop.name, 'w')

    seen_defaults = -999
    for i,line in enumerate(tmp_top_f.readlines()):
        if line.startswith("; nbfunc"):
            seen_defaults = i

        # nbfunc, comb-rule, gen-pairs, fudgeLQ, fudgeQQ
        if i == seen_defaults + 1:
            nbfunc,comb,genpairs,fudgelj,fudgeqq = tuple(line.split())

            if genpairs != 'yes' or fudgelj != '1':
                utils.printerr("WARNING: gen-pairs and fudgeLJ default parameters are incorrect for SepTop")
                utils.printerr("    Fixing gen-pairs and fudgeLJ for compatibility with SepTop")
                genpairs_st = 'yes'
                fudgelj_st = '0.5'
                line = f"{nbfunc:<16}{comb:<16}{genpairs_st:<16}{fudgelj_st:<13}{fudgeqq:<12}\n"                    
        
        tmp_modtop_f.write(line)
    tmp_modtop_f.close()
            
    new_structure = parmed.load_file(tmp_top.name, xyz=tmp_gro.name)
    new_structure.save(tmp_modtop.name, overwrite=True)

    # delete all temporary files
    os.unlink(tmp_gro.name)
    os.unlink(tmp_modtop.name)
    os.unlink(tmp_top.name)

    return new_structure

def unique_atom_types(pmd_structure, ligcode):
    ''' adapted from: Gaetan Calabro
    '''
    omm_system = pmd_structure.createSystem(nonbondedMethod=app.NoCutoff,
                                                constraints=None,
                                                removeCMMotion=False,
                                                rigidWater=False)
    topology = pmd_structure.topology
    positions = pmd_structure.positions
    def check_water(res):
        two_bonds = list(res.bonds())
        if len(two_bonds) == 2:
            waters = []
            for bond in two_bonds:
                elem0 = bond[0].element
                elem1 = bond[1].element
                if (elem0.atomic_number == 1 and elem1.atomic_number == 8) \
                            or (elem0.atomic_number == 8 and elem1.atomic_number == 1):
                    waters.append(True)
                if all(waters):
                    return True
        else:
            return False
        
    atom_types_dic = {}
    count_id = 0
    for c in topology.chains():
        for r in c.residues():
            for a in r.atoms():
                if r.name + a.name in atom_types_dic:
                    a.id = atom_types_dic[r.name + a.name]
                else:
                    if check_water(r):
                        if a.element.atomic_number == 1:
                            a.id = 'HW'
                        else:
                            a.id = 'OW'
                            atom_types_dic[r.name + a.name] = a.id
                    elif r.name == ligcode:
                        a.id = 'L' + str(count_id)
                        # else:
                        #     a.id = 'O' + str(count_id)
                        #     atom_types_dic[r.name + a.name] = a.id

                    count_id += 1

    new_system_structure = parmed.openmm.load_topology(topology,
                                                        system=omm_system,
                                                        xyz=positions)
    new_system_structure.positions = pmd_structure.positions
    new_system_structure.velocities = pmd_structure.velocities
    new_system_structure.box = pmd_structure.box
        
    return new_system_structure

def find_clashing_water(pmd_struct, lig_resname, distance):
    ''' Find waters that are sterically clashing with a ligand.
        
        Parameters:
        ===========
        pmd_struct : parmed.Structure
            The structure to analyze.
        lig_resname : str
            The up-to-three character residue name.
        distance : float
            The distance cutoff (in nanometers) for clash detection.
            
        Returns:
        ========
        water_resnums : Iterable[int]
            The residue numbers of waters that are clashing with the ligand.

        adapted from: Gaetano Calabro
        
    '''
    with tempfile.NamedTemporaryFile(suffix='.pdb') as tf:
        app.PDBFile.writeFile(pmd_struct.topology, pmd_struct.positions, open(tf.name, 'w'))
        traj = mdtraj.load(tf.name)
    top = traj.topology

    lig_atom_idxs = top.select(f'resname {lig_resname}')
    lig_res_idx = top.atom(lig_atom_idxs[1]).residue.index
    wat_atom_idxs = top.select('resname HOH and name O')
    wat_res_idxs = [top.atom(i).residue.index for i in wat_atom_idxs]
    potential_contacts = [(lig_res_idx, wat_res_idx) for wat_res_idx in wat_res_idxs]
    contacts = mdtraj.compute_contacts(traj,contacts=potential_contacts, scheme='closest',ignore_nonprotein=False)


    # Note that this is 0-indexed, while the parmed structure is 
    # 1-indexed, therefore we add 1 before returning
    clash_res_idx = [i[1]+1 for i in contacts[1][(contacts[0] < 0.1)[0,:]]]
    return clash_res_idx










