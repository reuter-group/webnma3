#! /usr/bin/env python

# DESCRIPTION: Command line interface for using webnma 
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)

import sys
import argparse
import shutil
from os import makedirs
from os.path import isfile, isdir, join, basename

from .config import DIR_DEFAULT_NAME, DIRS_S, DIRS_C, MODE_NM
from . import calc_nm, deformation, mode_vis, fluc_disp, overlap, correlation
from . import profile_alignment, similarity
from .utils.pdb import download_pdb, webnma_save
from .utils.modefiles import verify_modes
from .utils import mustang
from .utils.webnma_exceptions import Webnma_exception

# <postional argument> <help> <type>  # must be=3
modefile_arg = 'modefile; mode filename to load; readable_file\n'
protein_arg = 'protein_file; protein structure file(.pdb, .ent or .cif); readable_file\n'

# <optional flag> <optional argument> <help> <dest> <default> [type] [choices] # >=5
path_arg_optional = '-p; --path; output path, default: .; out_path; .; path\n '
    
nm_cmd = \
    'mode: Calculate normal modes\n' \
    + protein_arg \
    + path_arg_optional \
    + '-n; --number; mode number to compute, default: {0}; mode_num; {0}; int \n'.format(MODE_NM) \
    + '-o; --out; output filename, default: modes.txt; out_file; modes.txt'

eigen_cmd = \
    'eigen: Calculate average deformation energies and plot eigenvalues\n' \
    + modefile_arg + protein_arg + path_arg_optional

fluc_cmd = \
    'fluc: Calculate fluctuations and atomic displacement analysis\n' \
    + modefile_arg + protein_arg + path_arg_optional

vis_cmd = \
    'mode_vis: Calculate protein trajectories for specified mode numbers/range\n' \
    + modefile_arg + protein_arg + path_arg_optional + \
    '-l; --lower; lower mode number range, default:7; lower_mn; 7; int \n' \
    '-u; --upper; upper mode number range, default:12; upper_mn; 12; int \n' \
    '-a; --amplitude; amplitude for protein motion, default:12; ampl; 12; int; 4,8,12,16,20 \n' \
    '-t; --trajectory; trajectory/frame number for the motion movie, default:8 ; traj; 8 ; int; 4,8\n'

corr_cmd = \
    'corr: Perform correlation analysis \n' \
    + modefile_arg + protein_arg + path_arg_optional + \
    '-n; --negative; threshold for negative correlation sticks, '\
    'can be percent (suggested 95~100) or hard value(suggested -1 ~ -0.8), default: 99; '\
    'neg; 99; float \n' \
    '-t; --threshold; threshold for positive correlation sticks, '\
    'can be percent (suggested 95~100) or hard value(suggested 0.8 ~ 1), default: 99; '\
    'pos; 99; float \n'

# overlap_cmd = \
#     'overlap: Perform overlap analysis \n' \
#     + modefile_arg + protein_arg + \
#     'overlap_file; protein structure to overlap; readable_file\n' \
#     '-e; --export; export the transformed C-alpha structure of the overlap file, default: 1; ' \
#     'export; 1; int; 0,1 \n' \
#     + path_arg_optional


cmds = [nm_cmd, eigen_cmd, fluc_cmd, vis_cmd, corr_cmd] # overlap_cmd]


parser = argparse.ArgumentParser(
    description='WEBnma v3 normal mode calculation and analyses'
)
sub_parser = parser.add_subparsers(dest='cmd_name')


def readable_file(string):
    if isfile(string):
        return string
    else:
        msg = "%s is not an existing file." % string
        raise argparse.ArgumentTypeError(msg)

def path(string):
    if isdir(string):
        return string
    else:
        msg = "%s is not an existing path." % string
        raise argparse.ArgumentTypeError(msg)


type_mapping = {'int':int,
                'float':float,
                'readable_file': readable_file,
                'path': path,
}

    
def args_dict(args):
    arg_dict = dict(help=args[0], dest=args[1], default=args[2])
    if len(args) >= 4:
        arg_dict['type'] = type_mapping[args[3]]
    if len(args) >= 5:
        arg_dict['choices'] = \
        [type_mapping[args[3]](e) for e in args[4].split(',')]
    return arg_dict


for cmd in cmds:
    lines = cmd.strip().split('\n')
    cmd_name, cmd_help = lines[0].split(':')
    p = sub_parser.add_parser(cmd_name, help=cmd_help)
    for l in lines[1:]:
        args = [e.strip() for e in l.strip().split(';' + ' ')]
        if len(args) == 3:  # postional argument
            p.add_argument(args[0], help=args[1], type=type_mapping[args[2]])
        elif len(args) >=5:   # optional argument
            d = args_dict(args[2:])
            p.add_argument(*args[:2], **d)
        else:
            print('Unknown argument:')
            print(args)

            
single_cmd_h = 'Perform normal mode analyses for a single structure'
sa_parser = sub_parser.add_parser('sa', help=single_cmd_h)
sa_parser.add_argument('pdb', 
                       help='PDB ID or protein structure file(.pdb, .ent or .cif)')

sa_parser.add_argument('-m', '--modefile',
                       help='mode filename to use',
                       dest='sa_modefile',
                       type=readable_file,
)
sa_parser.add_argument('-c', '--corr',
                       help='perform correlation analysis',
                       dest='corr_flag', action='store_true')

sa_parser.add_argument('-o', '--overlap',
                       help='PDB ID or protein struture file for overlapping',
                       dest='overlap_pdb')

sa_parser.add_argument('-p', '--path', help='output path',
                       default='.', dest='out_path', type=path)

sa_parser.add_argument('-z', '--zip', help='zip the result files', dest='zip_name')

sa_parser.add_argument('-s', '--select', nargs='+', help='chain selection', dest='chains')

sa_parser.add_argument('-v', '--overlap_select', nargs='+',
                       help='chain selection for overlap structure',
                       dest='o_chains')


dl_cmd_h = 'Download protein structure files from PDBe; support chain selection'
dl_parser = sub_parser.add_parser('dl', help=dl_cmd_h)
dl_parser.add_argument('pdb_id', help='PDB ID')
dl_parser.add_argument('-f', '--format', help="default: pdb, or mmcif ",
                       default='pdb', dest='format')

dl_parser.add_argument('-s', '--select', nargs='+',
                       help='Chain selection, e.g., A B',
                       dest='chains')

dl_parser.add_argument('-p', '--path', help='output path', default='.',
                       dest='out_path', type=path)


# Comparative analysis
ca_cmd_h = 'Perform comparative analyses for multiple structures'
ca_parser = sub_parser.add_parser('ca', help=ca_cmd_h)
ca_parser.add_argument('pdbs',
                       nargs='+',  # can recieve >= 1 structures
                       help='Protein structure files(.pdb, .ent or .cif)')

ca_parser.add_argument('-s', '--similarity',
                       help='perform covariance similarity analysis',
                       dest='simi_flag', action='store_true')

ca_parser.add_argument('-a', '--alignment',
                       help='provide an alignment file',
                       dest='align_file')

ca_parser.add_argument('-p', '--path', help='output path',
                       default='.', dest='out_path', type=path)

ca_parser.add_argument('-z', '--zip', help='zip the result files', dest='zip_name')


def main():
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()
        
    cmd = parser.parse_args(sys.argv[1:])
    
    if cmd.cmd_name == 'mode':
        calc_nm.main(cmd.protein_file, cmd.out_path, cmd.out_file, cmd.mode_num)
    
    elif cmd.cmd_name == 'eigen':
        deformation.main(cmd.modefile, cmd.protein_file, cmd.out_path)
        
    elif cmd.cmd_name == 'fluc':
        fluc_disp.main(cmd.modefile, cmd.protein_file, cmd.out_path)
        
    elif cmd.cmd_name == 'mode_vis':
        if (cmd.lower_mn > 0) and (cmd.upper_mn > 0) \
             and (cmd.lower_mn <= cmd.upper_mn):
            mode_range = (cmd.lower_mn, cmd.upper_mn)
            mode_vis.main(cmd.protein_file, cmd.modefile, cmd.out_path,  
                          mode_range, cmd.ampl, cmd.traj)
        else:
            msg = "Usage Error: lower/upper range must be postive integers, and lower <= upper."
            print(msg)
                    
    elif cmd.cmd_name == 'corr':
        if (-1 <= cmd.neg <= -0.8 or 95 <= cmd.neg <= 100) and \
           (0.8 <= cmd.pos <= 1 or 95 <= cmd.pos <= 100):
            threshold = (cmd.neg, cmd.pos)
            correlation.main(cmd.modefile, cmd.protein_file, cmd.out_path,
                         threshold)
        else:
            msg = 'Usage Error: please input proper threshold.' 
            print(msg)
#            raise argparse.ArgumentTypeError(msg)
        
    elif cmd.cmd_name == 'overlap':
        overlap.main(cmd.modefile, cmd.protein_file, cmd.overlap_file,
                     cmd.out_path, cmd.export == 1)
        
    elif cmd.cmd_name == 'dl':
        if len(cmd.pdb_id) != 4:
            print("Invalid PDB ID. Please check your PDB ID.")
        else:
            pdb_file = download_pdb(cmd.pdb_id, cmd.out_path, cmd.format)
            if cmd.chains:
                new_pdb_file = pdb_file[:-4]+'_selected.pdb'
                webnma_save(pdb_file, new_pdb_file, chains=cmd.chains)
            
    elif cmd.cmd_name == 'ca':
        if len(cmd.pdbs) < 2:
            print("Usage Error: please input >1 structures for comparative analysis")
        else:
            # Creat the folder for storing job analysis results
            if cmd.zip_name is not None:
                job_dir = join(cmd.out_path, cmd.zip_name)
            else:
                job_dir = join(cmd.out_path, DIR_DEFAULT_NAME)
                
            makedirs(job_dir, exist_ok=True)

            # preprocess: use only C-alpha
            new_pdbs = []
            for pdb in cmd.pdbs:
                new_pdb_ca = join(job_dir, basename(pdb))
                if new_pdb_ca[-4:].lower() == '.cif':
                    new_pdb_ca = new_pdb_ca[:-4] + '.pdb' 
                try:
                    webnma_save(pdb, new_pdb_ca, c_alpha=True, max_size=10000)
                except Webnma_exception as e:
                    print(e.error_str)
                    sys.exit(e.exit_code)

                new_pdbs.append(new_pdb_ca)

            # Default analysis dir, i.e., profile 
            makedirs(join(job_dir, DIRS_C[0]), exist_ok=True)            
                        
            # Alignment
            if cmd.align_file is not None:
                align_fasta = join(job_dir, basename(cmd.align_file))
                shutil.copy(cmd.align_file, align_fasta)
            else:
                makedirs(join(job_dir, DIRS_C[0], 'visualization'), exist_ok=True)
                mustang.main(new_pdbs, job_dir, 'alignment')
                align_fasta = join(job_dir, "alignment.afasta")
                                 
            profile_alignment.main(align_fasta, join(job_dir, DIRS_C[0]))
            
            # Covariance similarity
            if cmd.simi_flag:
                makedirs(join(job_dir, DIRS_C[1]), exist_ok=True)                
                similarity.main(align_fasta, join(job_dir, DIRS_C[1]))

            if cmd.zip_name is not None:
                shutil.make_archive(job_dir, 'zip', job_dir)
                shutil.rmtree(job_dir)
                
    # default single analysis
    elif cmd.cmd_name == 'sa':
        # Creat the folder for storing job analysis results
        if cmd.zip_name is not None:
            job_dir = join(cmd.out_path, cmd.zip_name)
        else:
            job_dir = join(cmd.out_path, DIR_DEFAULT_NAME)          
        makedirs(job_dir, exist_ok=True)
        for d in DIRS_S:            
            makedirs(join(job_dir, d), exist_ok=True)
    
        if len(cmd.pdb) == 4:            
            pdb_file = download_pdb(cmd.pdb, tar_dir=job_dir)
        else:
            pdb_file = cmd.pdb
            
        ca_pdbfile = join(job_dir, 'modes_CA.pdb')

        # select C-alpha and chains
        try:
            webnma_save(pdb_file, ca_pdbfile, c_alpha=True, chains=cmd.chains, max_size=20000)
            if cmd.chains:
                # keep also a full PDB with selected chains for secondary structure calculation
                chain_selected_pdbfile = join(job_dir, basename(pdb_file)[:-4]+"_chain_selected.pdb")
                webnma_save(pdb_file, chain_selected_pdbfile, chains=cmd.chains)
        except Webnma_exception as e:
            print(e.error_str)
            sys.exit(e.exit_code)

        modefile = join(job_dir, 'modes.txt')  # use default modefile name
        if cmd.sa_modefile is None:
            calc_nm.main(ca_pdbfile, job_dir)
        else:            
            if verify_modes(cmd.sa_modefile, ca_pdbfile):
                shutil.copy(cmd.sa_modefile, modefile)
            else:
                print('WARNING: unsupported mode file format in %s; recompute the modes.' % cmd.sa_modefile)
                calc_nm.main(ca_pdbfile, job_dir)
        
        deformation.main(modefile, ca_pdbfile, join(job_dir, DIRS_S[0]))
        fluc_disp.main(modefile, chain_selected_pdbfile if cmd.chains else pdb_file, join(job_dir, DIRS_S[1]))
        mode_vis.main(ca_pdbfile, modefile, join(job_dir, DIRS_S[2]))

        if cmd.corr_flag:
            correlation.main(modefile, ca_pdbfile, join(job_dir, DIRS_S[3]))
    
        if cmd.overlap_pdb is not None:
            if len(cmd.overlap_pdb) == 4:
                o_pdb = download_pdb(cmd.overlap_pdb, join(job_dir, DIRS_S[4]))
            else:
                o_pdb = cmd.overlap_pdb
                
            o_pdb_new =  o_pdb[:-4]+"_selected.pdb"
            try:
                webnma_save(o_pdb, o_pdb_new, chains=cmd.o_chains, c_alpha=True)
            except Webnma_exception as e:
                print(e.error_str)
                sys.exit(e.exit_code)    
            overlap.main(modefile, ca_pdbfile, o_pdb_new, join(job_dir, DIRS_S[4]))

        if cmd.zip_name is not None:
            shutil.make_archive(job_dir, 'zip', job_dir)
            shutil.rmtree(job_dir)
            

if __name__ == '__main__':
    main()
    # print(parser.parse_args(sys.argv[1:]))
