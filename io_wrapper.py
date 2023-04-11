import tempfile

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
                    if len(resid_ct) != 0 and not keep_resid(resid_ct[0]) and not line.startswith(";"):
                        continue

                if line.startswith("[ molecules ]"):
                    seen_molecules = True
                    
                outfile.write(line)
    return

