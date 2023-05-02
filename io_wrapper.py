import tempfile
import subprocess

def rm_molecule(top, mod_top, keep_resids):
    ''' remove molecules listed in keep_resids (HOH, SOL, ions) from
        [ molecules ] section of .top file
        This function preserves the FF parameters for HOH/SOL and 
        ions but removes the molecule name and count from the
        [ molecules ] section
        
        
        Parameters:
        ===========
        top:          str    input .top file name with all molecules
                             in [ molecules ] section of .top file
        mod_top:      str    output .top file name without unwanted 
                             resids in [ molecules ] section of .top 
                             file
        keep_resids:  list   list of residue ids to keep
        
        
        Returns:
        ========
        None
        
    '''
    with open(top, 'r') as infile:
        with open(mod_top, 'w') as outfile:
            seen_molecules = False
            for line in infile:
                if seen_molecules:
                    def keep_resid(resid):
                        for r in keep_resids:
                            if resid == r:
                                return True
                        return False
                    
                    resid_ct = line.split()
                    if len(resid_ct) != 0 and not keep_resid(resid_ct[0]) and not line.startswith(";") \
                            or len(resid_ct) == 0:
                        continue

                if line.startswith("[ molecules ]"):
                    seen_molecules = True
                    
                outfile.write(line)
    return

def is_a_subset(sub_file, super_file, disregard_n_lines=2, cutoff=0.02):
    ''' sub_file is a subset of super_file2
    '''
    ct = 0
    with open(super_file, 'r') as super_f:
        with open(sub_file, 'r') as sub_f:
            super_lines = []
            for l in super_f.readlines():
                if l.strip() != "":
                    super_lines.append(l.rstrip()[:])
            sub_lines = []
            for l in sub_f.readlines():
                if l.strip() != "":
                    # :44 goes up to end of coordinates before velocities
                    # :20 goes up to end of residue/atom information
                    sub_lines.append(l.rstrip()[:44])

            if len(super_lines) != len(sub_lines):
                return False

            for i,line in enumerate(super_lines):
                if i < disregard_n_lines:
                    continue
                if not line.startswith(sub_lines[i]):
                   #print("first diff at", i+1, line, sub_lines[i])
                   return False

    return True


def is_approx_same(sub_file, super_file, disregard_n_lines=2, cutoff=0.003):
    ''' sub_file is a subset of super_file
    '''
    ct = 0
    with open(super_file, 'r') as super_f:
        with open(sub_file, 'r') as sub_f:
            super_lines = []
            for l in super_f.readlines():
                if l.strip() != "":
                    super_lines.append(l.rstrip())
            sub_lines = []
            for l in sub_f.readlines():
                if l.strip() != "":
                    # :44 goes up to end of coordinates before velocities
                    # :20 goes up to end of residue/atom information
                    sub_lines.append(l.rstrip()[:44])

            if len(super_lines) != len(sub_lines):
                print(f"diff number of lines thus diff numbers of atoms")
                return False

            for i,line in enumerate(super_lines):

                sub_line = sub_lines[i]

                if i < disregard_n_lines:
                    continue

                if i == 1: # atomct line:
                    if not line == sub_line:
                        print(f"diff number of atoms: {line} :: {sub_line}")
                        return False
                    else:
                        continue

                if i + 1 == len(super_lines):
                    if not line == sub_line:
                        print(f"diff last line: {line} :: {sub_line}")
                        return False
                    else:
                        continue

                supres,supatm,supaid,supx,supy,supz = tuple([line[0:8].strip(),
                                                            line[8:15].strip(),
                                                            line[15:20].strip(),
                                                            line[20:28].strip(),
                                                            line[28:36].strip(),
                                                            line[36:44].strip()])
                subres,subatm,subaid,subx,suby,subz  = tuple([sub_line[0:8].strip(),
                                                            sub_line[8:15].strip(),
                                                            sub_line[15:20].strip(),
                                                            sub_line[20:28].strip(),
                                                            sub_line[28:36].strip(),
                                                            sub_line[36:44].strip()])

                supx = float(supx); supy = float(supy); supz = float(supz)
                subx = float(subx); suby = float(suby); subz = float(subz)

                resid_atm_aid_same =  supatm ==  subatm #supres == subres and supatm == subatm and supaid == subaid

                x_close = subx < supx + cutoff and subx > supx - cutoff
                y_close = suby < supy + cutoff and suby > supy - cutoff
                z_close = subz < supz + cutoff and subz > supz - cutoff

                # if "LIG" in supres:
                #     print(subx, , supx - cutoff, supx + cutoff)

                if resid_atm_aid_same and x_close and y_close and z_close:
                    pass # line is approx same

                else:
                    print(f"diff @ line #{i+1}; {line} :: {sub_line}")
                    print(supatm, subatm)
                    print(resid_atm_aid_same, x_close, y_close, z_close)
                    return False

    return True



def diff(file1, file2):
    ''' returns the number of lines different between files 1 and 2
    '''
    result = subprocess.run(
        [f'diff --suppress-common-lines --side-by-side {file1} {file2} | wc -l'],
        shell=True,
        capture_output=True,
    )
    return int(result.stdout.strip())


def count_lines(file):
    result = subprocess.run(
        [f'cat {file} | wc -l'],
        shell=True,
        capture_output=True,
    )
    return int(result.stdout.strip())

def assert_atomtypes(top, mod_top):
    top_atomtypes = []
    modtop_atomtypes = []

    with open(top, 'r') as topf:
        seen_atomtypes = False
        
        for line in topf:
            if seen_atomtypes and len(line.strip()) == 0:
                break

            if line.startswith(';'):
                continue

            if seen_atomtypes:
                name,num,mass,_,_,_,_ = tuple(line.split())
                top_atomtypes.append((name,num,mass))

            if line.startswith("[ atomtypes ]"):
                seen_atomtypes = True


    with open(mod_top, 'r') as mod_topf:
        seen_atomtypes = False

        for line in mod_topf:
            if seen_atomtypes and len(line.strip()) == 0:
                break

            if line.startswith(';'):
                continue

            if seen_atomtypes:
                name,num,mass,_,_,_,_ = tuple(line.split())
                modtop_atomtypes.append((name,num,mass))

            if line.startswith("[ atomtypes ]"):
                seen_atomtypes = True

    if len(top_atomtypes) != len(modtop_atomtypes):
        return False, top_atomtypes, modtop_atomtypes

    diff_top_atomtypes = []
    diff_modtop_atomtypes = []
    for ta, mta in zip(top_atomtypes, modtop_atomtypes):
        if ta != mta:
            diff_top_atomtypes.append(ta)
            diff_modtop_atomtypes.append(mta)

    if len(diff_top_atomtypes) != 0:
        return False, diff_top_atomtypes, diff_modtop_atomtypes

    return True, [], []