# DESCRIPTION: Python wrapper for using mustang to align PDBs
#   e.g mustang -i 2vl0_A.pdb 3eam_A.pdb 3rhw_A.pdb 3tlu_A.pdb -o webnma -F fasta
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)

import subprocess
import os, sys
from os.path import join
import matplotlib.pyplot as plt
import seaborn as sb

from .webnma_exceptions import ALIGNMENT_VISUALIZATION_FAIL, MUSTANG_FAIL
# from .testing import timing


# @timing 
def main(pdbs, tar_dir=".", identifier='webnma3'):
    '''
    when the alignment is done successfully, the following will be generated:
    1. <identifier>.afasta : alignment file 
    2. <identifier>.html : alignment in HTML
    3. <identifier>.rms_rot : RMSD matrix
    4. <identifier>.pdb : superposition of all structures
    5. profile/visualiztion/str00, str01, ... strxx : alignmented structures
    6. rmsd_clustered.txt : clustered RMSD matrix data file 
    '''
    old_dir = os.path.abspath(os.curdir)

    # use relative path so the output HTML won't expose the absolute path
    os.chdir(tar_dir)
    pdbs_rel = [os.path.basename(p) for p in pdbs]
    
    pdbs_str = '-i ' + ' '.join(pdbs_rel)
    alignment_format = '-F fasta -F html'
    rmsd = '-r ON'
    output_identifier = '-o ' + identifier
    cmd = ' '.join(['mustang', pdbs_str, alignment_format, rmsd, output_identifier])
    
    cmd += ' | grep "Running" >> %s.html' % identifier  # add running time to the file

    # Step 1: alignment
    status1, _, err1 = launch_cmd(cmd)

    error_file = "mustang_error.log"
    if status1 != 0 :  # crucial error 
        with open(error_file, 'w') as f:
            f.write(err1.decode())
        sys.exit(MUSTANG_FAIL().exit_code())
        
    # Step2: calculate the cluster matrix from the raw RMSD matrix
    rmsd_cluster(identifier+".rms_rot", pdbs_rel)
    
    # Step3: split the identifier.pdb file by 'TER*' (and repeat this len(pdbs)-2 times)
    # for alignment visualization
    cmd2 = "csplit -s -f profile/visualization/str  %s.pdb  /TER*/+1 '{%s}'"\
        % (identifier, len(pdbs)-2)

    # Step3+ : add END after each chain so Pymol can process them as independent structures
    # NOTE: sed works differently on Mac OS VS Linux,
    # on Mac cmd line is: sed 's/^TER.*$/&\'$'\nEND/' %s.pdb
    cmd3 = "sed 's/^TER.*$/&\\'$'\nEND/' %s.pdb > %s_END_sep.pdb" % (identifier, identifier)
    cmd4 = 'mv %s_END_sep.pdb %s.pdb' % (identifier, identifier)
    
    status2, _, err2 = launch_cmd(cmd2 + ';' + cmd3 + ';' + cmd4)
    
    if status2 != 0 :   # not crucial, but warning  ??
        raise ALIGNMENT_VISUALIZATION_FAIL

    os.chdir(old_dir)
    

# Extract the RMSD matrix from mustang-generated .rms_rot file,
# then use seaborn.clustermap to sort/cluster this matrix,
# and restore the result to `rmsd_clustered.txt` for front-end RMSD heatmap visualization
def rmsd_cluster(filename, pdbs) -> None:
    mat_data = []
    pdbs = [p[:-4] if p[-4:] in ['.ent','.pdb','.cif'] else p for p in pdbs]
    
    with open(filename) as f:
        content =  f.readlines()
        mat_string = [line[4:].split() for line in content[3:3+len(pdbs)]]
        mat_data = [[float(s) if s != '---' else .0 for s in ss]
                         for ss in mat_string]
            
    ax = sb.clustermap(mat_data)
    plt.gcf().clear()
        
    reordered_pdbs = [pdbs[i] for i in ax.dendrogram_row.reordered_ind]
    with open("rmsd_clustered.txt",'w') as f:
        f.write('\t'.join(reordered_pdbs) + '\n')
        f.writelines(['\t'.join([str(d) for d in l])+'\n' for l in ax.data2d.values])


    
def launch_cmd(cmd: str):
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         shell=True,
                         stderr=subprocess.PIPE,
                         bufsize=1)

    stdout, stderr = p.communicate()
    p.wait()
    return (p.returncode, stdout, stderr)


if __name__ == '__main__':
    import sys
    from os.path import join
    if len(sys.argv) > 2:
        main(sys.argv[1:]) # do not pass identifier
    else:
        INPUT = join('webnma_api','tests', 'data_profile_alignment', 'input')
        OUTPUT = INPUT[:-5] + 'output'
        pdbs = ['2vl0_A.pdb', '3eam_A.pdb', '3rhw_A.pdb', '3tlu_A.pdb']
        pdbs_full = [join(INPUT, p) for p in pdbs]
        main(pdbs_full, OUTPUT)
