# DESCRIPTION: Analyses for correlation matrix
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)

from os.path import basename, join
from typing import List, Tuple
import numpy as np
import json 

from .utils.modefiles import load_modes
from .utils.pdb import calc_dist
from .utils.webnma_plot import plot_and_record
from .config import DIS_MIN

# some configurations for correlation sticks visualization
COLORS = ['red', 'blue']
ATOMS_USED = "CA"


def calc_corre(modefile):
    es, modes, rs = load_modes(modefile, load_res=True)
    es, modes = es[6:], modes[:, 6:]
    atom_num = len(rs)

    es_mat = np.diag(es**(-2))

    covariance = modes @ es_mat @ modes.transpose()
    
    corre = np.zeros((atom_num, atom_num))
    for i in range(atom_num):
        for j in range(i+1):
            corre[i,j] = sum([covariance[i*3, j*3],
                            covariance[i*3 + 1, j*3 + 1],
                            covariance[i*3 + 2, j*3 + 2]])
    # standardize
    corre_std = np.zeros((atom_num, atom_num))                           
    for i in range(atom_num):
        for j in range(i+1):
            corre_std[i,j] = corre[i,j] / (corre[i,i] * corre[j,j])**0.5
            corre_std[j,i] = corre_std[i,j]

    return corre_std



def pymol_group(res_pairs, g_name, color, pdb_name):    
    res_set = set(np.array(res_pairs).ravel()) #  + " and chain ...)"
    list_atoms = ", ".join(["(id {})".format(r) for r in res_set])
    creates = "create %s, %s and name %s and (%s)" \
                               % (g_name, pdb_name, ATOMS_USED, list_atoms)

    bond_formatter = "bond ({0} in ({3} and id {1})), ({0} in ({3} and id {2}))"
    bonds = [bond_formatter.format(g_name, p[0], p[1], pdb_name) for p in res_pairs]
    
    show_bonds = [
        "show sticks, %s" % g_name,
        "set stick_radius, 0.2, %s" % g_name,
        "color %s, %s" % (color, g_name),
        "hide cartoon, %s" % g_name,
        "\n",
    ]

    return [creates] + bonds + show_bonds


def thr_format(t:float) -> str:
    if 0 < t <=1 :
        return '({:.3f} < value < 1)'.format(t)
    elif -1 <= t < 0:
        return '(-1 < value < {:.3f})'.format(t)
    else:
        return '({}% < value < 100%)'.format(t)
    

def pymol_script(res, pdb_path):
    """
    Function to write a .pml file (for PyMol) showing sticks for
    correlated and anti-correlated atoms
    res_pairs :: ([(int, int)], [(int, int)])
    """
    
    pdb_name = basename(pdb_path)[:-4]
    thr_upper, thr_lower = res[0][0], res[1][0]
            
    header = [
        "# Pymol script for visualising correlation networks on input protein from WEBnma",
        "# To visualise: first load modes_CA.pdb (can be downloaded from WEBnma website) in Pymol,",
        "# then run this script file within the same session.",
        "# Objects are created where pairs of residues that have a correlation value within a",
        "# threshold and a minimum distance of {}Å, ".format(DIS_MIN),
        "# are represented by sticks forming uninterrupted networks"
        "# if they cover a larger region within the matrix.",
        "# Red sticks represent selected positive correlations " + thr_format(thr_upper),
        "# Blue sticks represent selected negative correlations "+ thr_format(thr_lower),
        "\n"]

    body = [
        "load %s" % basename(pdb_path),
        "hide everything",
        "set cartoon_trace, 1",
        "show cartoon, name %s" % ATOMS_USED,
        "\n"]

    for i, g_res in enumerate(res):
        thr, pairs = g_res
        if len(pairs) > 0:
            g_name = "group{}_{:.3f}_{}".format(i, thr, pdb_name)
            group = pymol_group(pairs, g_name, COLORS[i], pdb_name)
            body += group
    
    return (header + body + ['zoom'])

        
def filter_res(corre_mat, dist_mat, thresholds) -> \
    Tuple[Tuple[float, List[Tuple[int, int]]], Tuple[float, List[Tuple[int, int]]]]:
    
    thr_lower, thr_upper = thresholds
    cor_res =[]
    neg_cor_res = []
    
    l = len(corre_mat)

    # hard threshold
    if (-1 <= thr_lower <= 0) & (0 <= thr_upper <=1):
        for i in range(l):
            for j in range(i):
                if dist_mat[i,j] > DIS_MIN:
                    if thr_upper <= corre_mat[i,j] <= 1:
                        cor_res.append((i+1,j+1))
                    elif -1 <= corre_mat[i,j] <= thr_lower:
                        neg_cor_res.append((i+1,j+1))
                        
    # percentile threshold
    else:
        pos=[]
        neg=[]
        for i in range(l):
            for j in range(i):
                if dist_mat[i,j] > DIS_MIN:
                    if corre_mat[i, j] > 0:
                        pos.append((corre_mat[i, j], (i+1, j+1)))
                    elif corre_mat[i,j] < 0:
                        neg.append((corre_mat[i, j], (i+1, j+1)))
        pos.sort(key= lambda e: e[0])
        neg.sort(key= lambda e: e[0])

        thr_upper = np.percentile([v[0] for v in pos], thr_upper)
        thr_lower = np.percentile([v[0] for v in neg], 100 - abs(thr_lower))
        
        cor_res = [e[1] for e in pos if e[0] > thr_upper]
        neg_cor_res = [e[1] for e in neg if e[0] < thr_lower]
        
    return ((thr_upper, cor_res), (thr_lower, neg_cor_res))


#####  Helpers for verification
# Load a correlation matrix from webnma2 .dat file 
def read_v2_cor_data(filename):
    with open(filename) as f:
        c = f.readlines()
        import numpy as np
        c_arr = np.array([[float(r) for r in l.split(' ')[1:]] for l in c[1:]])
        return c_arr

def compare_mat(m1, m2, atol=1e-8) -> bool:
    if m1.shape != m2.shape:
        print("They have different shapes.")
        return False
    if not np.array_equal(m1, m2):
        if not np.allclose(m1, m2, atol):
            print('They are not equal or close.')
            return False
        else:
            print('They are very close.')
            return True
    else:
        print('They are identical.')
        return True
######


def litemol_json(corre, dist_mat, tar_dir):
    for t in range(95,100):
        thresholds = (t,t)
        (_, pos_cor), (_, neg_cor) = filter_res(corre, dist_mat, thresholds)
        cor_dict = {'positive': pos_cor, 'negative': neg_cor}
        fname = "correlation_{}_{}.json".format(thresholds[0], thresholds[1])
        with open(join(tar_dir, fname), 'w') as f:
            json.dump(cor_dict, f)

        
        
# TODO: NEED to map array/list indexes to real residue indices ??
def main(modefile, pdbfile, tar_dir='.', thresholds=(99, 99)):
    corre = calc_corre(modefile)
    atom_num = len(corre)
    plot_and_record(
        (range(atom_num), 'Residue index (in all selected chains)'),
        (range(atom_num), 'Residue index (in all selected chains)'),
        title = 'Correlation Matrix',
        tar_dir = tar_dir,
        name = 'correlation.png',
        style = 'heatmap',
        z = corre,
        data_filename='correlation_matrix.txt'
    )
    
    dist_mat = calc_dist(pdbfile)
    
    # write correlated residue pairs to files for LiteMol visualization    
    litemol_json(corre, dist_mat, tar_dir)

    # write PyMol script
    cor_res = filter_res(corre, dist_mat, thresholds)
    fname = "correlation_{}_{}".format(thresholds[0], thresholds[1])
    fname = join(tar_dir, fname)
    with open(fname + '.pml', 'w') as f:
        lines = pymol_script(cor_res, pdbfile)
        f.writelines('\n'.join(lines))

        
    
if __name__ == '__main__':
    import sys
    main(*sys.argv[1:])

