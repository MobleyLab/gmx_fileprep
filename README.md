# gmx_fileprep
Scripts for local gromacs input file preparation:
* solvation
* openFF forcefield swapping

## Requires:
Please use gmxfileprep.yaml as your environment. Your environment should use parmed v4 or newer. This version of parmed adds atomtyping checks so atomtypes are not combined haphazardly. This ***should*** prevent atomtyping issues...

```
conda env create -f gmxfileprep.yaml
```

## Warnings
This program uses `parmed`, if you run into strange gromacs errors when running, the first thing I'd suggest to check is that all your atomtypes in the .top file are correct
