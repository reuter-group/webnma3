# DESCRIPTION: Analyses for overlap
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)

import numpy as np
import os
from .utils.modefiles import load_modes
from .utils.pdb import read_pdb, webnma_save, rewrite_pdb
from .utils.superimpose import calc_trans, apply_trans
from .utils.webnma_plot import plot_and_record, record
from .utils.webnma_exceptions import *

def mass_dot_product(p, q, mass):
    return sum(np.ravel(p * q * mass[:, np.newaxis]))


def mode_norm(m, mass):
    m_w = mass_dot_product(m, m, mass)
    return m / np.sqrt(m_w)


# return None if the two proteins are different
def overlap_analyse(modefile, pdb, overlap_pdb, export_trans=False):
    es, virb_modes = load_modes(modefile, vibr = True)
    PDB_ntuple = read_pdb(pdb, unit_nm=True)
    O_PDB_ntuple = read_pdb(overlap_pdb, unit_nm=True)

    atom_num = len(PDB_ntuple.residues)
    if atom_num == len(O_PDB_ntuple.residues):
        for r1, r2 in zip(PDB_ntuple.residues, O_PDB_ntuple.residues):
            f = r1 == r2
        if not f:
            return None
    else:
        return None
            
    mass = PDB_ntuple.weight**2
    P, Q = PDB_ntuple.ca_coords, O_PDB_ntuple.ca_coords

    # find the transformation for superposing Q on P
    trans, _ = calc_trans(P, Q, mass)
    Q = apply_trans(*trans, Q)

    # normalized difference
    diff = mode_norm(Q - P, mass)
        
    ps = []
    for m in virb_modes.transpose():
        m_norm = mode_norm(m.reshape(atom_num, 3), mass)  # normalized mode
        projection = mass_dot_product(m_norm, diff, mass)
        ps.append(projection ** 2)
        
    ps_acc = np.add.accumulate(ps)    
    return (ps, ps_acc, Q) if export_trans else (ps, ps_acc)


# pdb 
def main(modefile, pdb, overlap_pdb, tar_dir = '.', export = True):
    try:
        res = overlap_analyse(modefile, pdb, overlap_pdb, export)
    except Webnma_exception as exc:
        # this exc does not abort the whole job, but only this single analysis
        print("WARNING: overlap analysis failed. " + exc.error_str)      
        return None

    if res is None: 
        return None   # Abort the computation
    
    if export:
        ps, ps_acc, Q = res
        Q = [a * 10 for a in Q]  # transform back to angstrom
        temp = overlap_pdb + '_ca'
        webnma_save(overlap_pdb, temp, c_alpha=True)
        rewrite_pdb(temp, Q, os.path.join(tar_dir, 'overlap_transformed.pdb'))
        os.remove(temp)
    else:
        ps, ps_acc = res
        
    plot_and_record(
        (range(len(ps)), 'Mode index'),
        (ps, 'Squared overlap'),
        name = 'overlap.png',
        tar_dir = tar_dir,
        header="index\tvalue"
    )
    
    plot_and_record(
        (range(len(ps_acc)), 'Mode index'),
        (ps_acc, 'Cumulative squared overlap'),
        name = 'overlap_acc.png',
        tar_dir = tar_dir,
        header="index\tvalue"
    )

    record(
        range(1,len(ps_acc)+1),
        zip(ps, ps_acc),
        tar_dir = tar_dir,
        data_filename = 'overlap_all.txt',
        header="index\toverlap\tcumulative"
    )
    

# if __name__ == '__main__':  
#     # testing data
#     mf = 'test/data_overlap/input/1su4_modes.dat'
#     pdb_file = 'test/data_overlap/input/1su4.pdb'
#     overlap_file = 'test/data_overlap/input/overlap_file.pdb'
#     main(mf, pdb_file, overlap_file, 'test/data_overlap/output/')
    
#     # import sys
#     #main(*sys.argv[1:])

