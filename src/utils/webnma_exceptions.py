# DESCRIPTION: Customized sys.exit code for integrating into the job runners of
#   WEBnma3 web server's backend
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)

import sys

class FASTA_INVALID(Exception):
    def __init__(self, exc='Fasta file can not be parsed.'):
        super().__init__(exc)
        # print(exc)
        sys.exit(4)


class MUSTANG_FAIL(Exception): 
    def __init__(self, exc='MUSTANG alignment failed.'):
        self.exc = exc
        self.code = 6
        super().__init__(exc)

    def error_reason(self):
        return self.exc
        
    def exit_code(self):
        return self.code
        
        
class FASTA_FORMAT_ERROR(Exception):
    def __init__(self, exc='Fasta file format does not satisfy the requirements'):
        self.exc = exc
        self.code = 7
        super().__init__(exc)
        # sys.exit(7)        

    def error_reason(self):
        return self.exc
        
    def exit_code(self):
        return self.code
    
    

class ALIGNMENT_VISUALIZATION_FAIL(Exception): 
    def __init__(self, exc='alignment visualization failed due to file preprocessing problem.'):
        self.exc = exc
        self.code = 8
        super().__init__(exc)

    def error_reason(self):
        return self.exc
        
    def exit_code(self):
        return self.code


## Suggestion: rewrite the exceptions above as subclass of Webnma_exception

class Webnma_exception(Exception):
    '''
    Use sys.exit(exc.exit_code) to exit with the correct code
    when handling with this exception class
    '''
    def __init__(self, error_str, exit_code):
        self.error_str = error_str
        self.exit_code = exit_code
        super().__init__(error_str)


class PDBFILE_INVALID(Webnma_exception):
    def __init__(self, pdb, error_str=""):
        super().__init__("Invalid PDB file {}: ".format(pdb) + error_str, 2)
        # sys.exit(2)


class PDB_DOWNLOAD_FAIL(Webnma_exception): 
    def __init__(self, pdb, error_str=''):
        super().__init__("Download {} failed. ".format(pdb) + error_str, 3)
        print(self.error_str)
        sys.exit(self.exit_code) # this exception will abort a job in any case


class PDB_OVERSIZE(Webnma_exception): 
    def __init__(self, limit_size, actual_size):
        super().__init__('Oversized PDB with {} residues. Max.:{}'.format(actual_size, limit_size), 5)
        # sys.exit(5)


class CONVERGENCE_ERROR(Webnma_exception):
    def __init__(self, error_str='Diagonalization computation does not converge'):
        super().__init__(error_str, 9)


class MODE_INVALID(Webnma_exception): 
    def __init__(self, error_str='Invalid modes: not exact 6 zeros are returned.'):
        super().__init__(error_str, 10)


class PDB_BOND_INVALID(Webnma_exception): 
    def __init__(self, r1, r2, dis):
        self.error_str = "Invalid bond: the input PDB has a non-standard distance between the Calpha atoms " +\
                        "of {} and {} of {:.6f} angstroms. Please resubmit with a minimized structure.".format(
                        r1, r2, dis)
        super().__init__(self.error_str, 11)


class PDB_COORDINATE_INVALID(Webnma_exception): 
    def __init__(self, invalid_count, total_count):
        self.error_str = "Invalid CA atom coordinates: {} out of all {} pairs of CA atoms have identical coordinates. ".format(invalid_count, total_count) +\
            " Remove those duplicated coordinates to run the calculation."
        super().__init__(self.error_str, 12)