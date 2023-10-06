# DESCRIPTION: Analyses for fluctuations and atomic displacement
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)

import numpy as np
from os.path import join

from .utils.modefiles import load_modes, rawmode_to_vibrational
from .utils.pdb import mass_protein, calc_ss
from .utils.webnma_plot import plot_and_record, record
from .config import Units_k_B, TMP, TAR_SS_H, TAR_SS_B, MODE_NM


def calc_disp(modefile, mode_num):
    mode_num = mode_num -1  
    es, modes, rs = load_modes(modefile, load_res=True)
    weight = mass_protein(rs, full=True)
    ca_num = len(rs)

    mode_vib = rawmode_to_vibrational(modes[:, mode_num], es[mode_num], weight)
    mode_vib = np.array(mode_vib)
    s = sum(mode_vib**2)
    disp_norm = [100 * sum(a_coord**2) / s for a_coord in mode_vib.reshape(ca_num,3)]
    return disp_norm
    


def calc_fluc(modefile:str, mode_nums=[]):
    es, modes, rs = load_modes(modefile, load_res=True)
    
    ca_num = len(rs)
    mass = mass_protein(rs, full=True, weighted=False)

    if len(mode_nums) > 0 and max(mode_nums) <= len(es) :
        mode_nums = [n-1 for n in mode_nums]  
        es_n = [es[n] for n in mode_nums]    
        modes_n = [modes[:,n] for n in mode_nums]
    else:
        # use all non-trivial modes 
        es_n = es[6:]
        modes_n = np.array(modes[:,6:]).transpose()
    
    # e**2 below to actually have eigenvalues (squares of the frequencies)
    f = lambda coord, e: sum(coord**2)/ e**2
    flucs1 = [[f(a, es_n[i]) for a in m.reshape(ca_num,3)]
                              for i, m in enumerate(modes_n)]
    flucs1_tr = np.array(flucs1).transpose()
    
    s = Units_k_B * TMP / (2. * np.pi) **2 
    flucs = [s * sum(f)/mass[i] for i, f in enumerate(flucs1_tr)]

    normalized_flucs = flucs_norm(np.array(flucs))
    print("\n****************************************************************************")
    print("NB! calc_fluc() with updated normalization, without squaring [2023 update]")
    print(f"Sanity check: sum of flucs for all residues is... {sum(normalized_flucs)}")
    print("****************************************************************************\n")
    return normalized_flucs


# final transformation above, normalizing values in [0,100]
def flucs_norm(arr):
    s = sum(arr)
    return 100 * arr / s


def calc_ss_spans(pdbfile):
    sec_strs= calc_ss(pdbfile)
    ids = list(range(1, len(sec_strs)+1))   
    ss_dict = {TAR_SS_B[0]:[], TAR_SS_H[0]: []}
    for tar_ss in [TAR_SS_B, TAR_SS_H]:
        id_ss = filter(lambda p: p[1] in tar_ss, zip(ids, sec_strs))
        id_ss = [i for (i,_) in list(id_ss)]
        if len(id_ss) > 0:
            id0 = id_ss[0]
            count = 1
            for i in id_ss[1:]:
                if (id0 + count) == i :
                    count += 1
                else:
                    ss_dict[tar_ss[0]].append((id0, id0+count-1))
                    id0 = i
                    count = 1
            if id0 != id_ss[-1]:
                ss_dict[tar_ss[0]].append((id0, id_ss[-1]))
        
    return ss_dict


def main(modefile, pdbfile, tar_dir='.'):
    # fluctuation analysis for all modes
    mode_nums=range(7, MODE_NM+7)

    # calculation secondary structure
    try:
        ss_dict = calc_ss_spans(pdbfile)
        style = 'fluc_ss'
    except:
        ss_dict = None
        style = 'default'
    
    fs = calc_fluc(modefile, mode_nums)
    plot_and_record(
        (range(len(fs)), 'Residue index (in all chains)'),
        (fs,'Fluctuation'),
        'Normalized fluctuations for modes from %d to %d (NB! 2023 update)' \
                        % (mode_nums[0], mode_nums[-1]),
        tar_dir = tar_dir,
        name = 'fluctuations.png',
        style = style,
        ss_dict = ss_dict,
        header="index\tfluctuation",
    )
    
    # displacement analysis for mode 7 to 12
    all_dis =[]
    for mode_num in range(7,13):
        dis = calc_disp(modefile, mode_num)
        all_dis.append(dis)
        plot_and_record(
            (range(len(dis)), 'Residue index (in all chains)'),
            (dis,'Displacement'),
            'Normalized squared atomic displacements for Mode %d' % mode_num,
            tar_dir = tar_dir,
            name = 'disp_mode_%d.png' % mode_num,
            style = style,
            ss_dict = ss_dict,
            header="index\tdisplacement",
        )                

    record(range(1, len(all_dis[0])+1),
           np.array(all_dis).transpose(),
           tar_dir,
           'displacements.txt',
           style=style,
           header="\t".join(["index"] + ["mode{}".format(m) for m in range(7,13)]),           
    )

    
if __name__ == '__main__':
    import sys
    main(sys.argv[1], pdbfile='../data/pdbs/1su4.pdb', tar_dir=sys.argv[2])
    # modefile = '../data/modes/1su4_modes.dat'
    # main(modefile)  
   
