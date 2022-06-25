---
sidebar_label: System
title: core.System
---

## System Objects

```python
class System()
```

Template object to contain System initialization settings

**Attributes**:

- `smiles` _str_ - SMILES string of the polymer (use * as connecting point)
  
- `density` _float_ - Density of the system
  
- `builder` _Builder_ - Builder (One of `EMC` or `PSP`)
  
- `natoms_total` _int_ - Total number of atoms in the system, one of
  this attribute and `nchains_total` has to be
  provided but not both (providing both will result
  in an error); default: `None`
  
- `nchains_total` _int_ - Total number of polymer chains in the system, one
  of this attribute and `natoms_total` has to be
  provided but not both (providing both will result
  in an error); default: `None`
  
- `natoms_per_chain` _int_ - Number of atoms per polymer chain, one of this
  attribute, `mw_per_chain`, and `ru_per_chain`
  has to be provided but not more than 1
  (providing more than 1 will result in an
  error); default: `None`
  
- `mw_per_chain` _int_ - Molecular weight of the polymer, one of this
  attribute, `natoms_per_chain`, and `ru_per_chain`
  has to be provided but not more than 1
  (providing more than 1 will result in an
  error); default: `None`
  
- `ru_per_chain` _int_ - Number of repeating unit per polymer chain, one of
  this attribute, `natoms_per_chain`, and
  `mw_per_chain` has to be provided but not more
  than 1 (providing more than 1 will result in an
  error); default: `None`
  
- `end_cap_smiles` _str_ - SMILES string of the end-cap unit for polymers
  ; default: `"*C"` (hint: put `"*[H]"` for end
  capping with -H)
  
- `data_fname` _str_ - File name of the output data file, which will be
  read in by LAMMPS
  [read_data](https://docs.lammps.org/read_data.html)
  command; default: `"data.lmps"`

### write\_data

```python
def write_data(output_dir: str = '.', cleanup: bool = True) -> None
```

Method to make LAMMPS data file (which contains coordinates and force
field parameters)

**Arguments**:

- `output_dir` _str_ - Directory for the generated LAMMPS data file
  ; default: `"."`
  
- `cleanup` _bool_ - Whether to clean up files other than the LAMMPS data
  file PSP generated
  

**Returns**:

  None

## SolventSystem Objects

```python
class SolventSystem(System)
```

Template object to contain System with solvents initialization settings

**Attributes**:

- `smiles` _str_ - SMILES string of the polymer (use * as connecting point)
  
- `solvent_smiles` _str_ - SMILES string of the solvent
  
- `ru_nsolvent_ratio` _float_ - The ratio of total number of repeating units
  in the system and total number of solvent
  molecules
  
- `density` _float_ - Density of the system
  
- `builder` _Builder_ - Builder (One of `EMC` or `PSP`)
  
- `natoms_total` _int_ - Total number of atoms in the system, one of
  this attribute and `nchains_total` has to be
  provided but not both (providing both will result
  in an error); default: `None`
  
- `nchains_total` _int_ - Total number of polymer chains in the system, one
  of this attribute and `natoms_total` has to be
  provided but not both (providing both will result
  in an error); default: `None`
  
- `natoms_per_chain` _int_ - Number of atoms per polymer chain, one of this
  attribute, `mw_per_chain`, and `ru_per_chain`
  has to be provided but not more than 1
  (providing more than 1 will result in an
  error); default: `None`
  
- `mw_per_chain` _int_ - Molecular weight of the polymer, one of this
  attribute, `natoms_per_chain`, and `ru_per_chain`
  has to be provided but not more than 1
  (providing more than 1 will result in an
  error); default: `None`
  
- `ru_per_chain` _int_ - Number of repeating unit per polymer chain, one of
  this attribute, `natoms_per_chain`, and
  `mw_per_chain` has to be provided but not more
  than 1 (providing more than 1 will result in an
  error); default: `None`
  
- `end_cap_smiles` _str_ - SMILES string of the end-cap unit for polymers
  ; default: `"*C"` (hint: put `"*[H]"` for end
  capping with -H)
  
- `data_fname` _str_ - File name of the output data file, which will be
  read in by LAMMPS
  [read_data](https://docs.lammps.org/read_data.html)
  command; default: `"data.lmps"`

### write\_data

```python
def write_data(output_dir: str = '.', cleanup: bool = True) -> None
```

Method to make LAMMPS data file (which contains coordinates and force
field parameters)

**Arguments**:

- `output_dir` _str_ - Directory for the generated LAMMPS data file
  ; default: `"."`
  
- `cleanup` _bool_ - Whether to clean up files other than the LAMMPS data
  file PSP generated
  

**Returns**:

  None

