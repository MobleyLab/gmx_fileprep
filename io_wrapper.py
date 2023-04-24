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

def is_a_subset(sub_file, super_file):
    ''' sub_file is a subset of super_file2
    '''
    ct = 0
    with open(super_file, 'r') as super_f:
        with open(sub_file, 'r') as sub_f:
            super_lines = []
            for l in super_f.readlines():
                if l.strip() != "":
                    super_lines.append(l.strip())
            sub_lines = []
            for l in sub_f.readlines():
                if l.strip() != "":
                    sub_lines.append(l.strip())

            if len(super_lines) != len(sub_lines):
                return False

            for i,line in enumerate(super_lines):
                if sub_lines[i] not in line:
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

