# DESCRIPTION: Analyses for covariance similarity, in comparative analyses
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)

import numpy as np
from os.path import basename, join
import matplotlib.pyplot as plt
import seaborn as sb
from scipy.linalg import eigh

from .utils.modefiles import load_modes
from .utils.superimpose import calc_trans
from .utils.pdb import read_pdb
from .utils.alignment_utils import parse_fasta, make_mask
# from .utils.testing import timing


DEFFRAC = .75
POOLFRAC = .95
def compute_bc(a_values, a_modes, b_values, b_modes):
    a_cutoff = get_index_for_coverage(a_values, POOLFRAC) # i.e. truncmodes1 
    a = spectral_assembly(a_values[:a_cutoff], a_modes.transpose()[:a_cutoff])
    a = divide_by_trace(a)

    b_cutoff = get_index_for_coverage(b_values, POOLFRAC) # truncmodes2
    b = spectral_assembly(b_values[:b_cutoff], b_modes.transpose()[:b_cutoff])
    b = divide_by_trace(b)

    q = (a + b) / 2
    q_values, q_modes = eigh(q)
     
    q_values = q_values[::-1]
    q_vectors = q_modes.transpose()[::-1]
    
    q_cutoff = get_index_for_coverage(q_values, DEFFRAC)
    q_values = q_values[:q_cutoff]
    q_vectors = q_vectors[:q_cutoff]

    a_contr = subspace_project(q_vectors, a)
    b_contr = subspace_project(q_vectors, b)
    q_contr = np.diag(q_values)

    db = bhattacharyya_flex_term(a_contr, b_contr, q_contr)

    # i.e. compmodes
    s = max(get_index_for_coverage(a_values, DEFFRAC),
            get_index_for_coverage(b_values, DEFFRAC))
    
    if db < 0:
        ndb = - np.sqrt(abs(db) / s)
    else:
        ndb = np.sqrt(db / s)

    return np.exp(-ndb)


def compute_rmsip(modes1, modes2, cutoff=10) -> float:
    """
    Root Mean Squared Inner Product of the vector sets
    """  
    vectors1 = modes1[:, 0:cutoff].transpose()
    vectors2 = modes2[:, 0:cutoff].transpose()
    
    squared_inner_product = 0
    for vector1 in vectors1:
        for vector2 in vectors2:
            inner_product = np.dot(vector1, vector2)
            squared_inner_product += inner_product * inner_product
    msip = squared_inner_product/len(vectors1)
        
    return np.sqrt(msip)


# modes:: 2D np.array, keeping the original shape from the modefile
# rot_tensor :: 3*3 np.array
def rotate_mode(modes, rot_tensor):
    rot_tensor = rot_tensor.transpose()
    rot_func = lambda m: m @ rot_tensor
    rot_modes = [np.apply_along_axis(rot_func,1, m.reshape(-1,3))
                 for m in modes.transpose()]
    res =  np.array(rot_modes).reshape(len(rot_modes), -1).transpose()
    return res 


# @timing
def compute_similarity(ca_coords, modes, values, similarity_measure='rmsip'):
    num = len(ca_coords)
    similarity_matrix = np.zeros((num, num))
    np.fill_diagonal(similarity_matrix, 1)
    
    for i in range(num):
        s1_coord = ca_coords[i]
        s1_modes = modes[i]
        
        for j in range(i+1,num):
            s2_coord = ca_coords[j]
            s2_modes = modes[j]
            
            (rot_tensor, _), _ = calc_trans(s1_coord, s2_coord)
            rotated_s2_modes = rotate_mode(s2_modes, rot_tensor)
            
            if similarity_measure == 'rmsip':
                similarity_matrix[i,j] = compute_rmsip(s1_modes, rotated_s2_modes)
            else:
                similarity_matrix[i,j] = compute_bc(values[i], s1_modes,
                                                    values[j], rotated_s2_modes) 
            similarity_matrix[j,i] = similarity_matrix[i,j]

    return similarity_matrix


##### Helpers from Webnma API V2
def bhattacharyya_flex_term(a, b, q):
    d_q = np.linalg.eigvalsh(q)
    d_a = np.linalg.eigvalsh(a)
    d_b = np.linalg.eigvalsh(b)

    #i.e. 0.5 * ( np.sum(np.log(d_q)) - .5 * ( np.sum(np.log(d_a)) + np.sum(np.log(d_b)) ) )
    return 0.5 * (sum(np.log( d_q / np.sqrt(d_a * d_b))) )


def subspace_project(basis, matrix):
    """
    Expresses the matrix in the (sub) space defined by the basis    
    basis: k vectors of size n
    matrix n*n matrix, k<=n    
    return: k*k matrix U'AU, where U has the basis as columns and A=matrix
    """
    u = np.array(basis.transpose())
    return u.transpose() @ matrix @ u


def divide_by_trace(matrix):
    if len(matrix.shape) == 1:
        return matrix/sum(matrix)
    else: # len(matrix.shape) == 2:
        return matrix/matrix.trace()

    
def get_index_for_coverage(values, coverage):
    if coverage < 1:
        cut_sum = sum(values) * coverage
        cum_var = values.cumsum()
        for i,cum in enumerate(cum_var):
            if cum > cut_sum:
                return i+1
    else:
        return len(values)
    

def spectral_assembly(values, vectors):
    """
    constructs a symmetric matrix from its eigenvectors and eigenvalues
    """            
    p = np.array(vectors)
    v = np.array(np.diag(values))
    
    return p.transpose() @ v @ p


def calc_cov_mat(values, modes, mask):
    """
    Adapted from API V2: _calculate_effective_covariance_matrix_()
    Calculates the effective covariance matrix for the aligned residues
    Effect of unaligned residues on aligned residues are included
    """
    # nontrivial_values, nontrivial_modes = remove_trivial(self.raw_values, self.raw_modes)
    values = values[6:]
    modes = modes[:, 6:]
    
    #obtain covariance matrix
    cov = spectral_assembly(1/values, modes.transpose())
    
    # #check if all or none are masked
    if sum(mask) == len(mask):
        return cov
    # elif len(mask[mask == True]) == 0:
    #     raise StructureComparisonException("No residues aligned")
        
    #partition hessian
    mask = np.array(mask).repeat(3)
    Cp, Cpq, Cqp, Cq = partition(cov, mask)        
    cov_eff = Cp - Cpq @ np.linalg.pinv(Cq) @ Cqp

    return cov_eff


def calc_cov_modes(modefile, mask):
    '''
    Adapted from API V2: _calculate_effective_covariance_modes()
    Calculate effective principal components
    '''
    values, modes = load_modes(modefile, freq=False)    
    covariance_matrix = calc_cov_mat(values, modes, mask)
    
    
    effective_pc_values, effective_modes = eigh(covariance_matrix)

    # TODO make this an independent function
    effective_pc_values = effective_pc_values[6:]
    effective_modes = effective_modes[:, 6:]

    # WHY using descending order of values and modes?
    reversed_modes = effective_modes.transpose()[::-1].transpose()
    return effective_pc_values[::-1], reversed_modes


def partition(matrix, mask):
    """
    partitins the matrix into four submatrices:
    the cells with indecies True in mask
    the cells i,j with mask[i] == True and mask[j] == False
    the cells i,j with mask[i] == False and mask[j] == True
    the cells with indecies False in mask
    """        
    n1, n2 = matrix.shape
    if n1 == n2 and len(mask) == n1:
        mp = matrix.compress(mask, 0).compress(mask, 1)
        mpq = matrix.compress(mask, 0).compress(~mask, 1)
        mqp = matrix.compress(~mask, 0).compress(mask, 1)
        mq = matrix.compress(~mask, 0).compress(~mask, 1)

        return np.array(mp), np.array(mpq), np.array(mqp), np.array(mq)
    else:
        return None

    
def sim_plot(pdb_names, mat, title="", tar_dir='.'):  
    mat = np.round(mat,2)
    ax = sb.clustermap(mat, linewidth=0.5, cmap="BuGn",
                    xticklabels=pdb_names,
                    yticklabels=pdb_names )
   
    plt.title(title)  
    fig_name = join(tar_dir, title+".png")
    plt.savefig(fig_name, format='png', dpi=200)
    plt.gcf().clear()

    
    data_file = join(tar_dir, title+".txt")    
    with open(data_file, 'w') as f:
        f.write('\t'.join(pdb_names) + "\n")
        f.writelines(['\t'.join([str(d) for d in l])+'\n' for l in mat])
    
    reordered_pdb_names =[pdb_names[i] for i in ax.dendrogram_row.reordered_ind]
    data_file2 = join(tar_dir, title+"_clustered.txt")    
    with open(data_file2, 'w') as f:
        f.write('\t'.join(reordered_pdb_names) + "\n")
        f.writelines(['\t'.join([str(d) for d in l])+'\n'  
                        for l in ax.data2d.values ])

# @timing
def main(alignment_file, tar_dir='.', mode_format='webnma3'):
    pdbs, seqs = parse_fasta(alignment_file, pdb_fullpath=True)

    # load protein structs
    pdb_structs = [read_pdb(p) for p in pdbs]

    # calculate aligment masks 
    masks = make_mask(seqs)

    # C-alpha coords after masking 
    masked_coords = [np.array([c for c, f in zip(s.ca_coords, m) if f]) \
                     for s,m in zip(pdb_structs, masks)]
    
    # calculate effective(aligned parts) covariance modes & values
    if mode_format.lower() == 'webnma2':
        modefiles = [pdb[:-4] + '_modes.dat' for pdb in pdbs]
    else:
        modefiles = [pdb[:-4] + '_modes_v3.txt' for pdb in pdbs]
        
    eff_pairs = [calc_cov_modes(md,m) for md, m in zip(modefiles, masks)]
    eff_values, eff_modes = list(zip(*eff_pairs))

    pdb_names = [basename(p)[:-4] for p in pdbs]
    
    for method in ['rmsip', 'bc']:
        sim_mat = compute_similarity(masked_coords, eff_modes, eff_values, method)
        sim_plot(pdb_names, sim_mat, title=method.upper(), tar_dir=tar_dir)

    
# if __name__ == '__main__':
#     import sys    

#     fasta = 'tests/data_similarity/input/mustang_alignment.fasta'
#     out = 'tests/data_similarity/output'
#     if len(sys.argv) < 2:
#         main(fasta, out)
#     else:
#         main(sys.argv[1], sys.argv[2])
