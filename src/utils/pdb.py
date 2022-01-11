# DESCRIPTION: Common functions for dealing with PDB files
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)

from collections import namedtuple
from os.path import join, basename
from typing import List
from urllib.request import urlretrieve
from urllib.error import URLError

import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from scipy.spatial.distance import cdist

from .webnma_exceptions import *
from .residue_mass import RES_MASS
from .residue_name import RES_NAME

DIG = 3 # digits kept for atom coordinates in pdb files

PDB_ntuple = namedtuple('Pdb_ntuple', 'ca_coords residues residues_full weight')

def is_ca(a):
    '''
    Check if an atom is a valid c-alpha atom
    Note: the full id of an atom is the tuple:
     (structure id, model id, chain id, `residue id`, atom name, altloc)
    where a `residue id` is:
       (hetero-flag, sequence identifier, insertion code)
    1. The hetero-flag: this is 'H_' plus the name of the hetero-residue
    (e.g. 'H_GLC' in the case of a glucose molecule), or 'W' in the case
    of a water molecule.
    2. The sequence identifier in the chain, e.g. 100
    3. The insertion code,
    For more: https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ

    Example: ('xxx', 0, 'A', (' ', 334, ' '), ('CA', ' '))
    '''    
    # 3 checks for C-alpha atoms: (1)correct atom name 
    # (2) empty hetero flag (3) valid residue name
    return a.get_id().upper() == 'CA' \
        and (a.get_full_id()[3][0]).split().__len__() == 0 \
        and a.get_parent().get_resname().upper() in RES_NAME.keys()


is_cif = lambda f: f[-4:].lower() == '.cif'


def set_parser(protein_file):
    '''
    Choose the correct parser according to the protein file's format
    '''
    if is_cif(protein_file):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser


# Support .pdb and .cif
def read_pdb(pdb, unit_nm = False, bl_check=False) -> namedtuple:
    '''
    Extract the coordinates of C-alpha atoms, residue names and the mass,
      pdb: PDB file 
      unit_nm: use nm as unit, otherwise angstrom as default
      bl_check: check if bond lengths are shorter than forcefield cutoff 0.277908 nm;
                this is mainly used as part of pre-processing for mode calculation
    Raise:
        Webnma_exception: PDBFILE_INVALID, PDB_BOND_INVALID
    '''

    parser = set_parser(pdb)
    try:
        models = parser.get_structure(pdb[:4], pdb)
    except Exception as e:
        raise PDBFILE_INVALID(pdb, " the structure can't be parsed")

    if len(models) > 1:
        print("WARNING: more than one models found, will only use the first model.")

    if len(models) == 0:
        raise PDBFILE_INVALID(pdb, "no valid models or chains")

    CAs = []
    Rs = []
    Rs_full = []
    for chain in models[0]:
        if chain.id.split().__len__() == 0:  # empty chain name
            # TODO: give warning for missing chain
            print("WARNING: empty chain found.".format(chain.id))
            chain_name = '?'
        else:
            chain_name = chain.id

        for res in chain:
            res_name = res.get_resname().upper()

            if res.id[0] != ' ':  # skip hetero residues or waters
                continue
            if not res_name in RES_NAME.keys():
                print("WARNING: invalid residue name found: {}, will be ignored.".format(res_name))
                continue            
            if not 'CA' in res:
                print("WARNING: invalid residue {}: no C-alpha atom".format(res_name))
                continue

            CA = res['CA'] # C-alpha atom
            CA_coord = CA.get_coord()/10 if unit_nm else CA.get_coord()
            res_fullname = chain_name + '.' + res_name + str(res.id[1]) # e.g, A.GLY123

            if bl_check and len(CAs) > 0: 
                dis = sum((CA_coord - CAs[-1])**2)**0.5 
                if unit_nm: 
                    if dis < 0.277908:
                        raise PDB_BOND_INVALID(Rs_full[-1], res_fullname, dis*10)                    
                elif dis < 2.77908:
                    raise PDB_BOND_INVALID(Rs_full[-1], res_fullname, dis)                    

            CAs.append(CA_coord)
            Rs.append(res_name)
            Rs_full.append(res_fullname)

    if len(CAs) == 0:
        raise PDBFILE_INVALID(pdb, "no valid atoms are selected.")
    else:
        w = mass_protein(Rs)
        return PDB_ntuple(
                    ca_coords = np.array(CAs, dtype=np.float64),
                    residues = np.array(Rs), # dtype=('U', 3)
                    residues_full = np.array(Rs_full), # dtype=('U', ?)
                    weight = np.array(w, dtype=np.float64)
                )    


def mass_protein(Rs, full=False, weighted=True):
    '''
    Retrieve the mass for each residue in 'Rs'
    Rs: list of 3-letter residue names 
    full: set True if residue name is in the format A.GLY123 (used modefile)
    weighted: set True if mass needs to be square-rooted

    TODO: use as only internal function  __mass_protein
    '''
    if full:
        Rs = [r[2:5] for r in Rs]

    Rs_mass = [RES_MASS[r.upper()] for r in Rs]

    if weighted:
        Rs_mass = [r**0.5 for r in Rs_mass]
    return Rs_mass
        


def calc_dist(pdbfile):
    '''
    Calculate the distance matrix of the CA atoms (unit: angstrom Å)
    '''
    CAs = read_pdb(pdbfile).ca_coords
    return cdist(CAs, CAs)
   

def calc_ss(pdbfile) -> List[str]:
    '''
    Calculate the secondary structure of the protein.
    Code Structure
    H 	Alpha helix (4-12)
    B 	Isolated beta-bridge residue
    E 	Strand
    G 	3-10 helix
    I 	Pi helix
    T 	Turn
    S 	Bend
    - 	None
    '''

    '''
    dssp_dict_from_pdb_file simply Popen 'mkdssp' and then deals with its output,
    src: http://biopython.org/DIST/docs/api/Bio.PDB.DSSP%27-pysrc.html#dssp_dict_from_pdb_file
    '''
    # keys :: [(chainid, res_id)], eg [('A', (' ', 12, ' ')), ...]
    ss_dict, keys = dssp_dict_from_pdb_file(pdbfile)  # this supports .cif also
    
    # make the resides' order consistent as it is in C-alpha file (i.e 'modes_CA.pdb')
    # to plot fluctuations with 2nd structure correctly
    parser = set_parser(pdbfile) 
    protein = parser.get_structure(pdbfile[:4], pdbfile)
    ss_list = []
    for a in protein.get_atoms():
        if is_ca(a):
            full_id = a.get_full_id()
            new_key = (full_id[2], full_id[3])
            if new_key in ss_dict:
                ss_list.append(ss_dict[new_key][1])
    return ss_list


def comp_pdb(pdb1, pdb2):
    pdb1 = read_pdb(pdb1)
    pdb2 = read_pdb(pdb2)
    for x,y in zip(pdb1.ca_coords, pdb2.ca_coords):
        test1 = np.array_equal(x,y)
    test2 = pdb1.residues == pdb2.residues
    return test1 and test2


# Only output .pdb, no .cif (but can use .cif as input also)
# because Bio.Python no longer maiatains MMCIFIO
def rewrite_pdb(pdb, CA_coords, new_name):
    parser = set_parser(pdb)
    protein = parser.get_structure(pdb[:4], pdb)

    i = 0
    for a in protein.get_atoms():
        if is_ca(a):
            try:
                a.set_coord([round(c, DIG) for c in CA_coords[i]])
                i = i + 1
            except IndexError:
                raise Exception('Unequal number of C-alpha atoms.')
  
    if i != len(CA_coords):
        raise Exception('Unequal number of C-alpha atoms.')

    io = PDBIO()
    io.set_structure(protein)
    io.save(new_name)


SERVER = "https://www.ebi.ac.uk/pdbe/entry-files/download/"
CIF_FORM_URL = SERVER + "{}.cif"
PDB_FORM_URL = SERVER + "pdb{}.ent"
def download_pdb(pdb_id, tar_dir='.', file_format='pdb'):
    '''
    Download pdb from PDBe(Protein Data Bank in Europe)
    Doc: https://www.ebi.ac.uk/pdbe/api/doc/
    If the download fails, the whole job fails (ie, sys.exit)
    '''
    pdb_id = pdb_id.lower()
    if file_format.lower() in ['cif', 'mmcif']:
        url = CIF_FORM_URL.format(pdb_id)
    else:
        url = PDB_FORM_URL.format(pdb_id)        
    dl_path = join(tar_dir, basename(url))
    try:
        urlretrieve(url, dl_path)
        print("Downloaded %s successfully." % pdb_id)
        return dl_path
    except URLError as e:
        raise PDB_DOWNLOAD_FAIL(pdb_id, e.reason)


def download_pdbs(pdb_ids:list, tar_dir='.', file_format='pdb'):
    return [download_pdb(i, tar_dir, file_format) for i in pdb_ids]

    
class Webnma_select(Select):
    def __init__(self, ca=False, chains=None):
        self.chains = "".join(chains) if chains is not None else chains
        self.ca = ca

    def accept_model(self, model):  # reserve only the first model
        return 1 if model.get_id() == 0 else 0 
        
    def accept_chain(self, chain):
        if self.chains is None:  # reserve all chains
            return 1        
        elif chain.get_id() in self.chains: # reserve selected chains
            return 1
        else:
            return 0

    def accept_atom(self, atom):
        if not self.ca: # reserve all atoms
            return 1
        elif is_ca(atom): # reserve only c-alpha
            return 1
        else:
            return 0

        
def webnma_save(source, target, chains=None, c_alpha=False, max_size=None):
    '''
    store the source (pdb structure) into target 
    with only selected chains or c-alpha atoms
    raise:
        Webnma_exception.PDBFILE_INVALID, PDB_OVERSIZE
    '''
    parser = set_parser(source)
    try:
        protein = parser.get_structure(source[:4], source)
    except Exception as e:
        raise PDBFILE_INVALID(source, " the structure can't be parsed")

    if max_size is not None:
        actual_size = len(list(protein.get_residues()))
        if actual_size > max_size:
            raise PDB_OVERSIZE(max_size, actual_size)
        
    io = PDBIO() # NOTE: Biopython only support PDB format output at this point
    io.set_structure(protein)
    io.save(target, Webnma_select(c_alpha, chains), preserve_atom_numbering=False)
    
    
# if __name__ == '__main__':
#     import sys
#     webnma_save('../../data/pdbs/1su4.pdb', '../test/1su4_ca.pdb')
