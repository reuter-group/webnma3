[![Python Package using Conda](https://github.com/reuter-group/webnma3/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/reuter-group/webnma3/actions/workflows/python-package-conda.yml)
# WEBnma3 API
This is the source code for calculating normal modes and performing further analyses used on [WEBnma web server](http://apps.cbu.uib.no/webnma3/).

- WEBnma v3 API is a refactoring work of WEBnma v2 API, which means:
  - v2 and v3 provide the same functions, but are implemented in different code structures. In particular, v2 uses Python2, and v3 uses Python3 with part of the code being parallelized to make use of modern multicore CPUs.  
  - the core function, normal mode calculation, in v2 is provided by [MMTK](https://github.com/khinsen/MMTK), while in v3 it is re-implemented with Python3.

- When using WEBnma v2 or v3, please cite: 
  > Tiwari SP, Fuglebakk E, Hollup SM, Skjærven L, Cragnolini T, Grindhaug SH, Tekle KM, Reuter N. WEBnm@ v2.0: Web server and services for comparing protein flexibility. BMC Bioinformatics. 2014; 15:427


# Install:

### For MacOS, Linux:

*[Optional]* If you are not using mamba instead of conda, we recommend running the following for much faster conflict-resolution and environment setup in conda.
  ```
  conda install -n base conda-libmamba-solver
  conda config --set solver libmamba
  ```

There are two simple setup alternatives. The first is to install with `conda`, using
```
conda install -c bioconda -c salilab -c ddx webnma
```

Alternatively, it is also easy to install from source: 

  1. Download the repository, _e.g._ using `git clone https://github.com/reuter-group/webnma3.git`.
  
  2. Navigate to the new directory that is created, and enter `conda env create --file environment.yml`.

     *NB!* If you are using an Apple M1 processor (for which some library dependencies are not yet available), you should probably run instead: `CONDA_SUBDIR=osx-64  conda env create --file environment.yml`.
  
  4. Still on that terminal window, activate the new environment and setup webnma:
  ```
  conda activate webnma3-env
  python setup.py install
  ```

Done! See the [usage instructions](#usage) below.

### For Windows users:

WEBnma does not have a `conda` build for Windows unfortunately, because some a key dependency (MUSTANG) do not have one. One way to get around this might be to use [WSL](https://docs.microsoft.com/en-us/windows/wsl/).  

You are welcome to contribute your solution here if any.


# Usage:
1. Use webnma via the web interface: [WEBnma web server](http://apps.cbu.uib.no/webnma3/)
2. Or use webnma on your own machine after installing webnma successfully. Type `webnma` from your command line and you will see the usage instruction:
```
$ webnma
usage: webnma [-h] {mode,eigen,fluc,mode_vis,corr,sa,dl,ca} ...

WEBnma v3 normal mode calculation and analyses

positional arguments:
  {mode,eigen,fluc,mode_vis,corr,sa,dl,ca}
    mode                Calculate normal modes
    eigen               Calculate average deformation energies and plot
                        eigenvalues
    fluc                Calculate fluctuations and atomic displacement analysis
    mode_vis            Calculate protein trajectories for specified mode
                        numbers/range
    corr                Perform correlation analysis
    sa                  Perform normal mode analyses for a single structure
    dl                  Download protein structure files from PDBe; support
                        chain selection
    ca                  Perform comparative analyses for multiple structures

optional arguments:
  -h, --help            show this help message and exit

```

To see how to use each sub-command/analysis you can input `webnma <sub-command> -h`, e.g `webnma mode -h` to see how to run a normal mode calculation.

More usage examples:
```
$ webnma dl 1su4     # download the file pdb1su4.ent from PDBe

$ webnma mode pdb1su4.ent  # calculate the (default 200) normal modes for this structure

$ webnma sa 1su4   # perform all the default analyses for the single structure 1su4
```


# Dependencies:
- Python packages:
  - numpy
  - scipy
  - numba
  - biopython
  - matplotlib
  - seaborn
- [MUSTANG](https://lcb.infotech.monash.edu/mustang/)

PDB files used in WEBnma 3 are downloaded from PDBe(Protein Data Bank in Europe).
