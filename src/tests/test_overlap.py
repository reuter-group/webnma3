from os import remove, listdir
from os.path import join, exists

from ..overlap import main
from ..utils.webnma_verify import verify


TESTDATA_DIR = 'data_overlap'
INPUT = join('tests', TESTDATA_DIR, 'input')
OUTPUT = join('tests', TESTDATA_DIR, 'output')
VERIFY = join('tests', TESTDATA_DIR, 'verify')


def test_main():
    input_pdb = join(INPUT, 'pdb3eam.ent')
    overlap_pdb = join(INPUT, 'pdb3tlu.ent')
    modefile = join(INPUT, 'modes.txt')

    main(modefile, input_pdb, overlap_pdb, OUTPUT)

    for f in ['overlap_transformed.pdb', 'overlap.txt', 'overlap_acc.txt']:
        output_file = join(OUTPUT, f)
        assert exists(output_file),"can not find output file."
        assert verify(output_file, join(VERIFY, f))
    
    # clean output files
    [remove(join(OUTPUT,f)) for f in listdir(OUTPUT) if '.gitkeep' not in f]
