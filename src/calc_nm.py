# DESCRIPTION: Core module for computing normal modes
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)

from os.path import join

import numpy as np
from numba import jit, prange
from scipy.linalg import eigh, LinAlgError
from scipy.sparse.linalg import eigsh, ArpackNoConvergence
from scipy.spatial.distance import cdist

from .utils.pdb import read_pdb
from .utils.modefiles import write_modefile
from .utils.webnma_exceptions import *
from .config import MODE_NM

# Note for diagonalization:
# eigh: 
#     Solve an ordinary or generalized eigenvalue problem for a complex
#     Hermitian or real symmetric matrix.
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html
#
# eigsh:
#     Find k eigenvalues and eigenvectors of the real symmetric square
#     matrix or complex hermitian matrix A.
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html


# for Konrad's model
CONST_A = 8.6e5
CONST_B = 2.39e5
CONST_C = 1.28e2
CONST_D = 0.4  # = 4 angstrom

ZERO = 1e-08 

@jit(nopython=True, cache=True, parallel=False, nogil=True)         
def konrad_force_cons(diss):
    n = diss.shape[0]  
    ks = np.zeros(diss.shape, dtype=np.float64)
    for i in range(n):
        for j in range(n):
            r = diss[i,j]
            if r < CONST_D:        
                ks[i,j] = CONST_A * r - CONST_B
            else:
                ks[i,j] = CONST_C * r**(-6)
    return ks

## set `parallel` True may cause multiprocessing problem used in 
## comparative analysis
@jit(nopython=True, cache=True) # parallel=False, nogil=True) 
def build_H(CAs, ks, diss, mass):
    n = CAs.shape[0]
    h = np.zeros((n,n,3,3), dtype=np.float64)
    for i in prange(n):
        for j in prange(n):
            if i != j:
                r = (CAs[i][0] - CAs[j][0], 
                     CAs[i][1] - CAs[j][1],
                     CAs[i][2] - CAs[j][2])
                k = ks[i,j]
                d_sq = diss[i,j]**2
                for x in range(3):
                    for y in range(3):
                        h[i,j,x,y] = - r[x] * r[y] * k / d_sq

                h[i,i] += h[i,j]
                h[i,j] = h[i,j] / mass[i] / mass[j]
        # i == j        
        h[i,i] = - h[i,i] / mass[i] / mass[i]
    return h 


def calc_modes(CAs, mass, n=MODE_NM):
    diss = cdist(CAs, CAs)  # distance matrix of all C-Alpha atoms
    zero_counts = np.sum(diss==0)
    if zero_counts > len(CAs):
        raise PDB_COORDINATE_INVALID((zero_counts-len(CAs)) /2, len(CAs))

    ks = konrad_force_cons(diss)
    H = build_H(CAs, ks, diss, mass) # build Hessian matrix
    h = np.concatenate(np.concatenate(H, axis=1), axis=1) # "reshape" H from (N,N,3,3) to (3N,3N)

    fetch_mode_num = min(n+6, 3*len(CAs))
    if n > 20:
        try: 
            e, v = eigh(h, eigvals=(0,fetch_mode_num-1), overwrite_a=True, check_finite=False )
        except LinAlgError:       
            raise CONVERGENCE_ERROR()  
    else: 
        try:            
            e, v = eigsh(h,k=fetch_mode_num,which='SA')
        except ArpackNoConvergence:       
            raise CONVERGENCE_ERROR()  
    
    if not np.allclose(e[:6], np.zeros((6,)), atol=ZERO) or e[6] <= ZERO:    
        raise MODE_INVALID()

    return e, v



def main(pdbfile, tar_dir='.', filename='modes.txt', mode_num=MODE_NM, exc_exit=True):
    try:
        PDB_ntuple = read_pdb(pdbfile, unit_nm=True, bl_check=True)   
        CAs = PDB_ntuple.ca_coords  
        e,v = calc_modes(CAs, PDB_ntuple.weight, mode_num)
    except Webnma_exception as exc:
        if exc_exit:
            print(exc.error_str)
            sys.exit(exc.exit_code)  # manually set the exit code when needed
        else:
            raise exc  # leave the exception unhandled

    modefile = join(tar_dir, filename)
    write_modefile(e,v,modefile, PDB_ntuple.residues_full)

    return modefile



if __name__ == '__main__':
    import sys
    main(sys.argv[1], sys.argv[2], sys.argv[1][:-4]+'_modes_v3.txt', sys.argv[3])
    
