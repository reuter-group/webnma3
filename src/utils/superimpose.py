
# DESCRIPTION: Superimpose two protein structures; adapted code from MMTK
# CONTRIBUTORS: 


import numpy as np

# atoms_coord :: 2D numpy.ndarray
# mass :: 1D numpy.ndarray


def center_of_mass(atoms_coord:np.ndarray, mass:np.ndarray) -> np.ndarray:
    return sum(atoms_coord * mass[:, np.newaxis])/sum(mass)


def rotation(v):
    rot = quaternion_rot()
    return np.dot(np.dot(rot, v), v)


def quaternion_rot():
    _rot = np.zeros((3,3,4,4))
    _rot[0,0, 0,0] =  1
    _rot[0,0, 1,1] =  1
    _rot[0,0, 2,2] = -1
    _rot[0,0, 3,3] = -1
    _rot[1,1, 0,0] =  1
    _rot[1,1, 1,1] = -1
    _rot[1,1, 2,2] =  1
    _rot[1,1, 3,3] = -1
    _rot[2,2, 0,0] =  1
    _rot[2,2, 1,1] = -1
    _rot[2,2, 2,2] = -1
    _rot[2,2, 3,3] =  1
    _rot[0,1, 1,2] =  2
    _rot[0,1, 0,3] = -2
    _rot[0,2, 0,2] =  2
    _rot[0,2, 1,3] =  2
    _rot[1,0, 0,3] =  2
    _rot[1,0, 1,2] =  2
    _rot[1,2, 0,1] = -2
    _rot[1,2, 2,3] =  2
    _rot[2,0, 0,2] = -2
    _rot[2,0, 1,3] =  2
    _rot[2,1, 0,1] =  2
    _rot[2,1, 2,3] =  2
    
    return _rot 


def calc_trans(P:np.ndarray, Q:np.ndarray, mass:np.ndarray = None):
    '''
    Calculate the transformation matrix that can be applied to Q to
    get the minimum (mass-weighted, if mass != None) RMS against P
    '''
    if mass is None:
        mass = np.full(len(P),1)
    
    N = len(mass)
    mass_w  = mass / sum(mass)
    
    ref_cms = center_of_mass(Q, mass)
    Q_ref = Q - ref_cms

    pos = sum(P * mass_w[:, np.newaxis])
    possq = sum(sum((P**2 + Q_ref**2) * mass_w[:, np.newaxis]))
    cross = (mass_w[:, np.newaxis] * P).transpose() @ Q_ref

    k = np.zeros((4, 4))
    k[0, 0] = -cross[0, 0]-cross[1, 1]-cross[2, 2]
    k[0, 1] = cross[1, 2]-cross[2, 1]
    k[0, 2] = cross[2, 0]-cross[0, 2]
    k[0, 3] = cross[0, 1]-cross[1, 0]
    k[1, 1] = -cross[0, 0]+cross[1, 1]+cross[2, 2]
    k[1, 2] = -cross[0, 1]-cross[1, 0]
    k[1, 3] = -cross[0, 2]-cross[2, 0]
    k[2, 2] = cross[0, 0]-cross[1, 1]+cross[2, 2]
    k[2, 3] = -cross[1, 2]-cross[2, 1]
    k[3, 3] = cross[0, 0]+cross[1, 1]-cross[2, 2]

    for i in range(1, 4):
        for j in range(i):
            k[i, j] = k[j, i]
    k = 2.*k
    for i in range(4):
        k[i, i] = k[i, i] + possq - sum(pos*pos)
    e, v = np.linalg.eigh(k)  # eigenvalues are in ascending order   
    v = v[:,0]     # minimum eigenvalue's vector
    if v[0] < 0: v = -v
    if e[0] <= 0.:
        rms = 0.
    else:
        rms = np.sqrt(e[0])
    rot_tensor = rotation(v)
    trans_vector = rot_tensor @ (-ref_cms) + pos
    return (rot_tensor, trans_vector), rms



def apply_trans(rot, trans, Q):
    return Q @ rot.transpose() + trans

