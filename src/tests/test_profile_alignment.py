import os
from os.path import join, exists

from ..profile_alignment import calc_nm_multip, calc_nm_wrapper, main
from ..utils.modefiles import comp_modefile
from ..utils.webnma_verify import verify
from ..utils.testing import timing

TESTDATA_DIR = 'data_profile_alignment'
INPUT = join('tests', TESTDATA_DIR, 'input')
OUTPUT = join('tests', TESTDATA_DIR, 'output')
VERIFY = join('tests', TESTDATA_DIR, 'verify')

@timing
def _calc_nm_multip(pdbs:list, tar_dir):    
    calc_nm_multip(pdbs, tar_dir)

@timing
def _calc_nm_seq(pdbs:list, tar_dir):
    for pdb in pdbs:
        calc_nm_wrapper((pdb, tar_dir))

        
def test_profile_alignment_calc_nm():
    '''
    report the time of mode calculation in sequencial VS using multiprocessing
    '''
    pdbs = [f for f in os.listdir(INPUT) if f[-4:] == '.pdb']
    pdbs_fullpath = [join(INPUT,p) for p in pdbs]
        
    for f in [_calc_nm_multip, _calc_nm_seq]:
        f(pdbs_fullpath, OUTPUT) #timing
        for pdb in pdbs:
            out = join(OUTPUT, pdb[:-4] + '_modes_v3.txt')
            verify_file = join(VERIFY, pdb[:-4]+'_modes.dat')
            assert comp_modefile(out, verify_file)

@timing 
def test_main():
    '''
    include 1.mode calculation 2.fluctuation 3. deformation energies 
    '''
    alignment_file = join(INPUT, 'mustang_alignment.fasta')
    fluc_png = join(OUTPUT, 'alignment_fluctuations.png')
    fluc_data_file = join(OUTPUT, 'alignment_fluctuations.txt')
    verify_file = join(VERIFY, 'fluctuations_use_v3_mode_DIG6.txt')  
    
    main(alignment_file, OUTPUT)

    assert exists(fluc_png), "fluctuation png does not exist."
    assert exists(fluc_data_file), "fluctuation data file does not exist."
    assert verify(fluc_data_file, verify_file), "fluctuations are not equal."

    def_png = join(OUTPUT, 'deformation_energies.png')
    def_data_file = join(OUTPUT, 'deformation_energies.txt')
    def_verify_file = join(VERIFY, 'mustang_alignment.fasta_deformation_energy.dat') 
    assert exists(def_png), "deformation profile png does not exist."
    assert exists(def_data_file), "deformation profile data does not exist."
    # the values of energy are large so have to use a loose verification ?
    assert verify(def_data_file, def_verify_file, rtol=1.5e-04), "deformation profile results are not equal."

    