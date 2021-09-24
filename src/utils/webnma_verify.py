# DESCRIPTION: helper functions for verifying outputs of webnma3 vs 2:
#   (line by line) numbers should be equal, and words should be identical
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)


import re
import numpy as np

DIG = 3

def is_number(s):
    return re.match(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?$", s) is not None

convert_num = lambda s: round(float(s), DIG)  # also works with 'nan'


def equal_nums(l1, l2, atol=1e-05, rtol=1e-05):
    if len(l1) != len(l2):
        return False

    l1 = [convert_num(s) for s in l1]
    l2 = [convert_num(s) for s in l2]
    
    return np.allclose(np.array(l1), np.array(l2), atol=atol, rtol=rtol, equal_nan=True)


# Check if two files are identical or very similar
# note: not suitable for compare modefiles from WEBnma v2 VS v3 because:
# modes can not be identical: they can have opposite directions; the first
# 6 modes are even meanless to compare; use utils.modefiles.comp_modefile instead
def verify(file1, file2, ignore_text=False, atol=1e-05, rtol=1e-05):
    if not file1[-4:] in ['.txt','.pdb', '.ent','.dat']:
        print('Skipped non-supported files: %s, %s' %(file1, file2))
        return True
    
    f1 = open(file1)
    f2 = open(file2)
    c1 = f1.readlines()
    c2 = f2.readlines()

    len1 = len(c1)
    len2 = len(c2)

    if len1 != len2:
        f1.close()
        f2.close()
        print('%s and %s have different numbers of lines: %d, %d' \
              %(file1, file2, len1, len2))
        return False
    
    diff_lines = []
    for i, (l1,l2) in enumerate(zip(c1,c2)):
        l1 = l1.strip()
        l2 = l2.strip()
        l1_ls = l1.split()
        if is_number(l1_ls[0]):
            l2_ls = l2.split()
            if not equal_nums(l1_ls, l2_ls, atol=atol, rtol=rtol):                
                diff_lines.append(i)
        else:
            if not ignore_text:
                if not l1 == l2:
                    diff_lines.append(i)
    f1.close()
    f2.close()
    
    if len(diff_lines) == 0:
        # print('{1} is consider the same as {0}'.format(file1, file2))
        return True
    else:
        print('{1} is different from {0}'.format(file2, file1))
        print('{:.2f}% lines are different. '.format(len(diff_lines)/len1*100))
        diff0 = diff_lines[0]
        print('First different line number: %d' % (diff0+1))
        print("In {}: Line {} \n {}".format(file1, diff0+1, c1[diff0]))
        print("In {}: Line {} \n {}".format(file2, diff0+1, c2[diff0]))
        return False        


if __name__ == '__main__':
    import sys
    verify(sys.argv[1], sys.argv[2])
