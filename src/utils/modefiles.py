# DESCRIPTION: Common functions for loading/writing WEBnma mode(and eigenvalue) files
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)

import numpy as np
from .pdb import mass_protein, read_pdb

HEADER = \
    "Normal modes text file from WEBnma\n" + \
    "Row 5: Eigenvalues(including first 6 trivial modes)\n" +\
    "Row 6: Sequence according to input PDB file (e.g. A.Gly123, " + \
        "where A is chain, Gly is Glycine and 123 is position in sequence)\n" + \
    "Row 7 to last row: Normal modes vectors (columns: modes index, rows:" + \
        "vectors X1,Y1,Z1, ..., Xn,Yn,Zn where n is the number of atoms)"


DELM = ' '  # delimiter used in modefiles 
DIG = 6 # digits kept for eigenvalues or modes in modefiles


def omega_to_lambda(w):
    """
    Switch omega(raw eigenvalues in modefiles) squared values to 
    lambda(frequencies)
    """
    return abs(w)**0.5 / (2 * np.pi)


T = 300
Units_k_B = 1.3806513e-23 * 6.0221367e23 / 1.0e3 
# Transform a raw mode (used in MMTK/WEBnma2) into a Vibrational mode
def rawmode_to_vibrational(mode_arr, freq, weights):
    amplitude = np.sqrt(2 * T * Units_k_B) / (2 * np.pi * freq)
    v_mode = [coord * amplitude / weights[i//3]  for i, coord in enumerate(mode_arr)]
    return np.array(v_mode)


# support WEBnma2 & 3 formats and any others that fulfill the following rules
'''
Checking rules:
1. the last 3*N lines are considered mode data (where N is C-alpha atom number)
   and one mode resides one column (delimited by white space)
2. eigenvalues (in one line) are somewhere above the mode data 
3. residues (in one line) have full information, i.e., <chain>.<RES><NUM>
'''
def verify_modes(filename, pdbfile) -> bool:
    pdb_info = read_pdb(pdbfile)
    ca_num = len(pdb_info.residues_full)
    
    with open(filename) as f:
        content = f.readlines()
        if len(content) < ca_num * 3:
            return False
    
        modes = content[-(ca_num*3):]
        mode_num = len(modes[0].strip().split(DELM))
        for m in modes[1:]:
            ls = m.strip().split(DELM)
            if len(ls) != mode_num:
                return False
            for e in ls:
                try:
                    e = float(e)
                except ValueError:
                    return False

        eigenvalues = None
        residues_full = None
        pdb_residues = DELM.join(pdb_info.residues_full).lower()
        for l in content[:-(ca_num * 3)]:       
            if l.lower().strip() == pdb_residues:
                residues_full = l  
            else:
                ls = l.strip().split(DELM)
                if len(ls) == mode_num:
                    for e in ls:
                        try:
                            e = float(e)
                        except ValueError:
                            return False
                    eigenvalues = ls
        return (eigenvalues is not None) and (residues_full is not None)

        

# Load all eigenvalues and modes **including** the first 6 trivial ones
# Modes are in 2D array keeping the same shape/position as they are in the modefile
# i.e., Mode 6 is modes[:,6]
# note: default load frequencies, rather than raw eigenvalues
# note: specify format 'webnma2' will keep only 3 digits
def load_modes(filename, freq=True, vibr = False, load_res=False, format='webnma3'):
    f = open(filename)
    content = f.readlines()  # including '\n' also
    content = [l.strip() for l in content]
    f.close()
    
    if format[:6] == 'webnma':
        eigenvalues = [float(x) for x in content[4].split(DELM)]
        modes = [[float(v) for v in row.split(DELM)] for row in content[6:]]
        if format == 'webnma2':
            eigenvalues = [round(e, DIG) for e in eigenvalues]
            modes = [[round(v, DIG) for v in row] for row in modes]
    else:
        raise Exception('Unsupported modefile format: ' + format)
    
    if freq:
        eigenvalues = [omega_to_lambda(e) for e in eigenvalues]
    
    residues_full = None 

    if vibr:
        if not freq:
            eigenvalues = [omega_to_lambda(e) for e in eigenvalues]
        residues_full = [r for r in content[5].split(DELM)]
        weights = mass_protein(residues_full, full=True)
        modes = np.array(modes)
        for i in range(6, len(eigenvalues)):
            modes[:,i] = rawmode_to_vibrational(modes[:,i], eigenvalues[i], weights)
    
    if load_res:
        if residues_full is None:
            residues_full = [r for r in content[5].split(DELM)]
        return np.array(eigenvalues), np.array(modes), np.array(residues_full)
    else:
        return np.array(eigenvalues), np.array(modes)

    
# check if the modes in two modefiles are equal:
# eigenvalues must be numerically equal, modes can be equal or have different directions
def comp_modefile(f1, f2, atol=1e-02):
    e1,m1 = load_modes(f1, False, format='webnma2')
    e2,m2 = load_modes(f2, False, format='webnma2')
    
    l = min(len(e1), len(e2))
    e1 = e1[:l]
    e2 = e2[:l]
    test1 = np.allclose(e1, e2, atol=atol)
    if test1:
        print('Eigenvalues are equal.')
    else:
        print('Eigenvalues are NOT equal.')
        return False

    l = min(m1.shape[1], m2.shape[1])
    m1 = m1[:, 6:l]
    m2 = m2[:, 6:l]
    
    test2 = np.allclose(m1,m2,atol=atol)
    if test2:
        print('Modes are equal.')
        return True
    else:
        test3 = np.allclose(abs(m1),abs(m2),atol=atol)
        if test3:
            print('Modes have equal scalars, but different directions.')
            return True
        else:
            print('Modes are NOT equal.')
            return False
    
        

def write_modefile(eigenvalues, modes, filename, residues):
    array2str = lambda arr: [str(a) for a in arr]
    es = np.round(eigenvalues, DIG)
    content = [ HEADER,
                DELM.join(array2str(es)),
                DELM.join(residues)]
    header = "\n".join(content)
    fmt = "%.{}f".format(DIG)
    np.savetxt(filename, modes, header=header, comments="", fmt=fmt, delimiter=DELM)

    
if __name__ == '__main__':
    f = 'test_modefile.txt'
    write_modefile(np.array([1,2,3]), np.array([[4],[5],[6]]), f, 'xxx')
    print(load_modes(f))
