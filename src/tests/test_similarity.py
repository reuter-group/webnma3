from os import listdir, remove
from os.path import join

from ..similarity import main
from ..utils.webnma_verify import verify
# from ..utils.testing import timing


TESTDATA_DIR = 'data_similarity'
INPUT = join('tests', TESTDATA_DIR, 'input')
OUTPUT = join('tests', TESTDATA_DIR, 'output')
VERIFY = join('tests', TESTDATA_DIR, 'verify')


def test_main():
    fasta = join(INPUT, 'mustang_alignment.fasta')
    main(fasta, OUTPUT, 'webnma2')

    for method in ['RMSIP', 'BC']:
        v3_file = join(OUTPUT, method + '.txt')
        v2_file = join(VERIFY, 'mustang_fasta_V2_{}.dat'.format(method))

        assert verify(v3_file,v2_file, ignore_text=True)

    [remove(join(OUTPUT,f)) for f in listdir(OUTPUT) if '.gitkeep' not in f]
