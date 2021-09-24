
from os import listdir, remove
from os.path import join, exists

from ..calc_nm import main
# from ..utils.webnma_verify import verify 
from ..utils.modefiles import comp_modefile

TESTDATA_DIR = 'data_calc_nm'
INPUT = join('tests', TESTDATA_DIR, 'input')
OUTPUT = join('tests', TESTDATA_DIR, 'output')
VERIFY = join('tests', TESTDATA_DIR, 'verify')

def test_main():
    for pdb in listdir(INPUT):
        pdb_fullpath = join(INPUT, pdb)
        if pdb[-4:] == '.pdb':
            mode_file = pdb[:-4] + '_modes_v3.txt'
            outfile = main(pdb_fullpath, OUTPUT, mode_file)
            verify_file = join(VERIFY, pdb[:-4] + '_modes.dat')
            if exists(verify_file):   # perfer webnma2 modefile
                assert comp_modefile(outfile, verify_file)
            elif exists(verify_file[:-4] + '_v3.txt'):  # otherwise use webnma3 modefile
                assert comp_modefile(outfile, verify_file[:-4] + '_v3.txt')

            remove(outfile) # clean up

