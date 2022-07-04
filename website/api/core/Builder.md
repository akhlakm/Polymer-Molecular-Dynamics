---
sidebar_label: Builder
title: core.Builder
---

## EMC Objects

```python
class EMC(Builder)
```

Object to perform system structure generation using
[EMC](http://montecarlo.sourceforge.net/): Enhanced Monte Carlo package.
This object should be used as input argument of `System` or `Lammps`
objects

**Attributes**:

- `force_field` _str_ - Force field, options are `"pcff"`, `"opls-aa"`,
  `"opls-ua"`, and `"trappe"`

## PSP Objects

```python
class PSP(Builder)
```

Object to perform system structure generation using
[PSP](https://github.com/Ramprasad-Group/PSP): Polymer Structure Predictor
package. This object should be used as input argument of `System` or
`Lammps` objects

**Attributes**:

- `force_field` _str_ - Force field, options are `"opls-lbcc"`,
  `"opls-cm1a"`, `"gaff2-gasteiger"`, and `"gaff2-am1bcc"`
  
- `packmol_nloop` _int_ - Maximum number of optimization loops of Packmol
  (PSP uses Packmol to pack molecules into a box); default: None
  
- `packmol_precision` _float_ - Packmol avoids atom overlaps by ensuring
  a 2.0 Angs atom distance, this parameter determines how close the
  solution must be to the desired distances to be considered correct
  ; default: None

