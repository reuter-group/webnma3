
from os.path import join

from ..correlation import calc_corre, read_v2_cor_data, compare_mat


TESTDATA_DIR = 'data_correlation'
INPUT = join('tests', TESTDATA_DIR, 'input')
# OUTPUT = join('tests', TESTDATA_DIR, 'output') # no output 
VERIFY = join('tests', TESTDATA_DIR, 'verify')

# test by comparing the correlation matrices of V3 with V2
def test_calc_corre():
    for pdb_id in ['1abi', '1su4']:
        # use V2 modefiles to get a much closer result to V2's correlation matrix 
        modefile = join(INPUT, pdb_id + '_modes.dat')
        corre = calc_corre(modefile) 
        v2_mat = read_v2_cor_data(join(VERIFY, pdb_id + '_correlation_matrix.dat'))
        assert compare_mat(corre, v2_mat)
