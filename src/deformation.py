# DESCRIPTION: Analyses for deformation and eigenvalue output
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)

from os.path import join
from scipy.spatial.distance import cdist
import numpy as np

from .utils.modefiles import load_modes
from .utils.pdb import read_pdb
from .utils.webnma_plot import plot_and_record, record
from .config import DEFORMATION_MODE_NUM, PLOT_EIGEN_NUM


## Deformation energy calculation constants
FACTOR = 46402
FC_LENGTH = 0.7
DISTANCE_CUTOFF = 1.2


def deformation_atom_sum(coords, vib_modes):
    energies = deformation(coords, vib_modes)
    return [sum(e) for e in energies]


def deformation(coords, vib_modes):
    '''
    coords.shape = (N, 3), where N is the number of residues/C-alpha atoms
    vib_modes.shape = (3N, DEFORMATION_MODE_NUM)  # raw shape from the modefile

    Return a list of normalized deformation energies of each atom for each mode
    i.e. np.array(energies_list).shape = (DEFORMATION_MODE_NUM, N) 
    '''
    N = len(coords)
    Rij_mat = cdist(coords, coords)  # pair distance matrix
    Rij_mat_sq = Rij_mat**2   

    modes_num = vib_modes.shape[1]
    vib_modes = vib_modes.transpose().reshape(modes_num, N, 3)

    energies_list = []
    for v in vib_modes:
        energy_mat = np.zeros((N, N))
        for i in range(N):
            for j in range(i):
                if Rij_mat[i][j] <= DISTANCE_CUTOFF:
                    k = np.exp( (0.01 - Rij_mat_sq[i][j]) / FC_LENGTH**2 )
                    l = np.dot(v[i] - v[j], coords[i] - coords[j])**2 
                    energy_mat[i][j] = k * l / Rij_mat_sq[i][j]
                    energy_mat[j][i] = energy_mat[i][j]

        v_sum = np.sum(v**2)
        energy = FACTOR * sum(energy_mat) / 2 / v_sum 
        energies_list.append(energy)

    return energies_list



def main(modefile, pdbfile, tar_dir=''):
    es, modes = load_modes(modefile, vibr=True)    
    PDB_ntuple = read_pdb(pdbfile, unit_nm=True) 

    des = deformation_atom_sum(PDB_ntuple.ca_coords, modes[:, 6:6+DEFORMATION_MODE_NUM])
    
    plot_defor_num = min(DEFORMATION_MODE_NUM, len(es)-6)
    plot_and_record(
        (range(6, plot_defor_num+6), "Mode index"),
        ( des, "Average deformation energy"),
        tar_dir = tar_dir,
        title='Average deformation energies for the %d lowest-frequency modes' % plot_defor_num,
        name='deformation.png',
        style='bar',
        header='index\tenergy', # \t separated headers 
    )
   
    plot_eigen_num = min(PLOT_EIGEN_NUM, len(es)-6)
    plot_and_record(
        (range(0, plot_eigen_num+6), "Mode index"),
        (es[:plot_eigen_num+6], "Eigenvalues"),
        tar_dir = tar_dir,
        title='Eigenvalues for the %d lowest-frequency modes' % plot_eigen_num,
        name='eigen.png',
        style='bar',
        header='index\teigenvalue' # \t separated headers 
    )
    
    
if __name__ == '__main__':
    import sys
    main(sys.argv[1], sys.argv[2])
