# gmx_fileprep
This repo contains scripts for local gromacs input file preparation, tailored for MobleyLab/SeparatedTopologies:
* solvation: see [prep_complex.py](prep_complex.py)
    * gmx options:
        * `gen-pairs = no`
        * `fudgeLJ = 0.5`
      
* openFF small molecule forcefield swapping
    * solvent swap ex: see [swap_solvent_ligand_ff.py](swap_solvent_ligand_ff.py)
    * complex swap ex: see [swap_complex_ligand_ff.py](swap_complex_ligand_ff.py)

## Requires:
Please use gmxfileprep.yaml as your environment. Your environment should use parmed v4 or newer. This version of parmed adds atomtyping checks so atomtypes are not combined haphazardly. This ***should*** prevent atomtyping issues...

```
conda env create -f gmxfileprep.yaml
```

## Warnings
This program uses `parmed`, if you run into strange gromacs errors when running, the first thing I'd suggest to check is that all your atomtypes in the .top file are correct
