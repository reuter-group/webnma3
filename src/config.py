
# DESCRIPTION: Default configurations for webnma3
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)

MODE_NM = 200 # default mode numbers to calculate

PLOT_EIGEN_NUM = 50  # number of eigenvalues to plot

DEFORMATION_MODE_NUM = 14 # default number of modes for calculating deformation energy

TMP = 300  # temperature

Units_k_B = 1.3806513e-23 * 6.0221367e23 / 1.0e3 

# target secondary structures (helix and β-sheet) for fluctuation plots
TAR_SS_B = ['E']   # β-sheet, min: 2 res 
TAR_SS_H = ['G','H','I']    # helixes, min: 3 res
TAR_SS = TAR_SS_B + TAR_SS_H

# Correlation
DIS_MIN = 8.   # 8 angstrom, minimum distance threshold

# result dir template
DIR_DEFAULT_NAME = 'WEBnma_analyses_results' 
# single PDB analyses
DIRS_S = ('deformation',
          'displacement',
          'visualization',
          'correlation',
          'overlap')
DIRS_C = ('profile', 'covariance') # comparative analyses

# number of available processes when using multiprocessing
PROCESS_NUM = 8