# DESCRIPTION: generate the trajectory files(.pdb format) for a structure
#   with a specified mode number
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)


from os.path import join
import numpy as np

from .utils.modefiles import load_modes, rawmode_to_vibrational
from .utils.pdb import read_pdb, rewrite_pdb


def vis(pdbfile, modefile, mode_range=(7,12), ampl=12, traj=8, tar_dir='', prefix=''):

    pdb_ntuple = read_pdb(pdbfile)
    es, modes = load_modes(modefile)
    
    if ampl*2 % traj != 0:
        raise Exception('Improper amplitude or trajectory number')

    step = int(ampl*2 / traj)
    coeffs = range(-ampl, ampl+1, step)
    
    atom_num = len(pdb_ntuple.residues)
    for i in range(mode_range[0]-1, mode_range[1]):  # switch to 0-indexed
        m = rawmode_to_vibrational(modes[:,i], es[i], pdb_ntuple.weight)

        m = m.reshape(atom_num,3)
        for coeff in coeffs:
            vector = m * coeff    
            output_coords = pdb_ntuple.ca_coords + vector * 10

            #file naming: prefix_modeNumber_ampl_traj_frameNumber.pdb            
            output_format = [str(i+1),   # switch back to 1-indexed
                            'ampl'+str(ampl),
                            'traj'+str(traj),
                            str((coeff+ampl)//step)]
            if prefix:
                output_format  = [prefix] + output_format
                
            output_name = join(tar_dir, '_'.join(output_format)+ '.pdb')
            
            rewrite_pdb(pdbfile, output_coords, output_name)


def pack_pdb(tar_dir,mode_range, ampl, traj):
    '''
    Pack all the trajectory pdbs(for the same mode number) into one
    with 'MODEL' numbers, so it can be visualizad in PyMol as a movie
    '''
    for n in range(mode_range[0]-1, mode_range[1]):
        pack_name = "{}_ampl{}_traj{}_movie.pdb".format(n+1, ampl, traj)
        output = join(tar_dir, pack_name)
        with open(output,'w') as pack_f:
            frame_num = list(range(traj+1)) + list(range(traj,-1,-1))
            for i in frame_num:
                pack_f.write('MODEL  {}\n'.format(i))
                traj_pdb = output[:-9] + str(i) + '.pdb'
                traj_f = open(traj_pdb)
                pack_f.write(traj_f.read())
                pack_f.write('ENDMDL \n')
                traj_f.close()
        

# mode_range is 1-indexed, both ends included
def main(pdb, mode, tar_dir, mode_range=(7,12), ampl=12, traj=8, prefix=''):
    vis(pdb, mode, mode_range, ampl, traj, tar_dir, prefix)
    pack_pdb(tar_dir, mode_range, ampl, traj)

    
if __name__ == '__main__':
    import sys
    main(sys.argv[1], sys.argv[2], sys.argv[3], (6,7))
