# DESCRIPTION: Analyses for profiling alignment, in comparative analyses
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)

import os
from os.path import join, dirname, exists
import numpy as np
from multiprocessing import Pool
from math import inf
from Bio import SeqUtils

from . import calc_nm, fluc_disp, deformation
from .utils.modefiles import load_modes
from .utils.webnma_exceptions import *
from .utils.webnma_plot import plot_and_record
from .utils.alignment_utils import parse_fasta
from .utils.pdb import read_pdb
from .config import PROCESS_NUM
# from .utils.testing import timing


def calc_nm_wrapper(pair):
    pdb, tar_dir = pair
    modefile_name = os.path.basename(pdb[:-4]) + '_modes_v3.txt'
    exc_code = 0
    try:
        calc_nm.main(pdb, tar_dir, modefile_name, mode_num=inf, exc_exit=False) #calculate all modes        
    except Webnma_exception as exc: 
        #! the child's exception can't be passed to or detected by parent process properly in Pool, 
        # so only pass the exit_code as a simple solution
        print("Mode calculation exception happened to %s :\n%s" % (pdb, exc.error_str))
        exc_code = exc.exit_code 

    return (os.path.join(tar_dir, modefile_name), exc_code)


# calculate modes for each pdb in parallel(use multiprocessing)
# @timing
def calc_nm_multip(pdbs:list, tar_dir='.'):
    ps = [(p, tar_dir) for p in pdbs]
    with Pool(processes=min(len(pdbs), PROCESS_NUM)) as pool:
        results = pool.map(calc_nm_wrapper, ps)
        for r, e in results:
            if e != 0:                 
                sys.exit(e) #! any mode calculation failure will abort the whole job
      
    return [r for (r, _) in results]


def deformation_profile(mode_pdb_pair):
    modefile, pdbfile = mode_pdb_pair
    _, vib_modes = load_modes(modefile, vibr=True)    
    PDB_ntuple = read_pdb(pdbfile, unit_nm=True)
    energies = deformation.deformation(PDB_ntuple.ca_coords, vib_modes[:, 6:])
    e_sum = sum(energies) * len(PDB_ntuple.ca_coords) / len(vib_modes[0][6:])
    return e_sum


def scatter(arr, seq)-> np.ndarray:
    i = 0
    new_arr = np.full(len(seq), np.nan)
    for j, a in enumerate(seq):
        if a == '-':
            continue
        else:
            new_arr[j] = arr[i]
            i += 1
    return new_arr

# @timing
def main(alignment_file, tar_dir='.'):
    pdbs, seqs = parse_fasta(alignment_file)
    work_dir = dirname(alignment_file)
    
    # fasta file validity (against pdb files) check
    # 1. all seqs have the same length
    # 2. pdb file names are correct
    # 3. sequences (in pdb and fasta) are identical
    valid = all([len(s) == len(seqs[0]) for s in seqs[1:]])

    if valid:   
        pdbs_fullpath = [] 
        for (pdb,seq) in zip(pdbs, seqs):
            pdb_fullpath = join(work_dir, pdb)
            if exists(pdb_fullpath):
                pdb_structure = read_pdb(pdb_fullpath)
                
                # convert 1-letter amino acid to 3-letter, "MetGluAlaAlaHis..."            
                seq_res = SeqUtils.seq3(list(filter(lambda x: x!= '-', seq.upper() )))
                
                if ''.join(pdb_structure.residues) != seq_res.upper():
                    valid = False
                    break
                else:
                    pdbs_fullpath.append(pdb_fullpath)
            else:
                valid = False
                break
    
    if not valid:
        sys.exit(FASTA_FORMAT_ERROR().exit_code())

    # compute all modes using multiprocessing
    # NOTE: modefiles are stored in the same dir of alignment
    modefiles = calc_nm_multip(pdbs_fullpath, dirname(alignment_file))
    
    with Pool(processes=min(len(pdbs), PROCESS_NUM)) as pool:
        # compute fluctuations for all non-trivial modes
        fs = pool.map(fluc_disp.calc_fluc, modefiles)
        # fs = [fluc_disp.calc_fluc(m) for m in modefiles]

        # Deformation energy profile
        def_es = pool.map(deformation_profile, zip(modefiles, pdbs_fullpath))
        # def_es = [deformation_profile(p) for p in zip(modefiles, pdbs_fullpath)]

    # scatter the fluctuation values according to alignment
    new_fs = [scatter(f, seq) for f,seq in zip(fs, seqs)]
    new_fs_tr = np.array(new_fs).transpose() # transpose for plotting
    plot_and_record(
        (range(len(seqs[0])), 'Alignment index (in all chains)'),
        (new_fs_tr,'Fluctuation'),
        'Normalized fluctuations for all modes',
        tar_dir = tar_dir,
        name = 'alignment_fluctuations.png',
        style = 'fluc_ss',
        ss_dict = {},
        pdbs = pdbs,
        header="\t".join(["index"]+[p[:-4]for p in pdbs])
    )

    new_def_es = [scatter(e, seq) for e, seq in zip(def_es, seqs)]
    new_def_es_tr = np.array(new_def_es).transpose() # transpose for plotting
    plot_and_record(
        (range(len(seqs[0])), 'alignment index'),
        (new_def_es_tr,'deformation energy'),
        'Normalized deformation energies averaged over all modes',
        tar_dir = tar_dir,
        name = 'deformation_energies.png',
        style = 'fluc_ss',
        ss_dict = {},
        pdbs = pdbs,
        header="\t".join(["index"]+[p[:-4]for p in pdbs])
    )


if __name__ == '__main__':
    import sys    
    # To run only this script: 
    # python -m src.profile_alignment [fasta_file]
    fasta = 'src/tests/data_profile_alignment/input/mustang_alignment.fasta'
    if len(sys.argv) < 2:
        main(fasta)
    else:
        main(sys.argv[1])
