# COBRA utils

[![PyPI version](https://badge.fury.io/py/cobra-utils.svg)](https://pypi.org/project/cobra-utils/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3470217.svg)](https://doi.org/10.5281/zenodo.3470217)

This is a series of utilities to use together with [COBRApy](https://github.com/opencobra/cobrapy).
The goal of this tool is to make easier the usage of COBRApy.

## Installation

To install this package from a Unix OS, run:

```
pip install cobra-utils
```

To install this package from Windows, run:
```
python -m pip install cobra-utils
```

In case you downloaded the source code, in a Unix OS run:

```
pip install -e .
```

or in Windows:
```
python -m pip install -e .
```

## Examples
* [Loading a model and retrieving reactions information](./notebooks/Ecoli_Rxn_Information.ipynb)
* [Reporter metabolites and pathways from differential expression in two strains of *E. coli*](./notebooks/Ecoli_Reporter_Metabolites_Pathways.ipynb)