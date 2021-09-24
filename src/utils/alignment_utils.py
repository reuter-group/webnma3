# DESCRIPTION: common function for processing FASTA file
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)

from os.path import join, dirname
from typing import Tuple
from Bio import SeqIO

from .webnma_exceptions import *

def parse_fasta(alignment_file, pdb_fullpath=False) -> Tuple[list, list]:
    pdbs = []
    seqs = []
    try:
        records = SeqIO.parse(alignment_file,'fasta')
        for r in records:
            pdb = r.id.split('|')[0]
            pdbs.append(pdb)
            seqs.append(r.seq) # r.seq::Bio.Seq.Seq, works as list-like
    except:        
        raise FASTA_INVALID("Fasta file can not be parsed.")
    
    if len(pdbs) == 0:
        raise FASTA_INVALID("Invalid fasta file.")
    
    if pdb_fullpath:
        pdbs = [join(dirname(alignment_file),pdb) for pdb in pdbs]
        
    return (pdbs, seqs)



def make_mask(seqs) -> list:
    '''
    seqs :: [seq1, seq2,...]
    for a set of aligned sequences, compute the conservation flag,
    i.e, a list of bools where a T represents that the residue is
    conserved across all the seqs
    '''
    conserved = lambda ls, check: sum([a != '-' for a in ls]) == check

    num_seq = len(seqs)
    flags = [conserved(col, num_seq) for col in zip(*seqs)]
    masks = [[f for (s,f) in zip(seq,flags) if s != '-' ] for seq in seqs]
    return masks
