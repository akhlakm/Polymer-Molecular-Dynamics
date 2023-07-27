import os
import re
import shutil
from io import TextIOWrapper
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import pyemc

from pmd.util import HiddenPrints, Pmdlogging, build_dir

PSP_FORCE_FIELD_OPTIONS = ('opls-lbcc', 'opls-cm1a', 'gaff2-gasteiger',
                           'gaff2-am1bcc')
EMC_FORCE_FIELD_OPTIONS = ('pcff', 'opls-aa', 'opls-ua', 'trappe')

EMC_EXTS = ('esh', 'data', 'in', 'params', 'vmd', 'emc.gz', 'pdb.gz', 'psf.gz',
            'cmap')
EMC_COEFF_EXCLUSIONS = ('bb', 'ba', 'mbt', 'ebt', 'at', 'aat', 'bb13', 'aa')


class Builder:

    def __init__(self, force_field: str,
                 force_field_options: Tuple[str]) -> None:
        self._force_field = force_field
        self._validate_force_field(force_field_options)

    def __repr__(self) -> str:
        return type(self).__name__

    def _validate_force_field(self, options):
        if self._force_field not in options:
            raise ValueError(f'Invalid {self} force_field, valid options are '
                             f'{", ".join(options)}')

    def write_data(self, output_dir: str, smiles: str, density: float,
                   natoms_total: int, length: int, nchains: int,
                   data_fname: str, cleanup: bool) -> None:
        raise NotImplementedError

    def write_solvent_data(self, output_dir: str, smiles: str,
                           solvent_smiles: str, density: float,
                           natoms_total: int, length: int, nsolvents: int,
                           nchains: int, data_fname, cleanup) -> None:
        raise NotImplementedError

    def write_functional_form(self, f: TextIOWrapper) -> None:
        raise NotImplementedError


class EMC(Builder):
    '''Object to perform system structure generation using
    [EMC](http://montecarlo.sourceforge.net/): Enhanced Monte Carlo package.
    This object should be used as input argument of `System` or `Lammps`
    objects

    Attributes:
        force_field (str): Force field, options are `"pcff"`, `"opls-aa"`,
            `"opls-ua"`, and `"trappe"`
    '''

    def __init__(self, force_field: str) -> None:
        super().__init__(force_field, EMC_FORCE_FIELD_OPTIONS)

    @staticmethod
    def _remove_brackets_around_asterisks(smiles: str) -> str:
        smiles = smiles.replace('[*]', '*')
        return smiles

    def _run_emc(self, tmp_file_prefilx: str, output_dir: str, data_fname: str,
                 cleanup: bool):
        Pmdlogging.info('Launching EMC...')
        try:
            previous_dir = os.getcwd()
            os.chdir(output_dir)
            pyemc.setup(f'{tmp_file_prefilx}.esh')
            pyemc.build('build.emc')

            # store parameters from .params file
            params = {}
            key = None
            with open(f'{tmp_file_prefilx}.params', 'r') as lines:
                for line in lines:
                    parts = line.split()
                    first_word = parts[0] if len(parts) > 0 else ''
                    # put an empty list when reaching a Coeffs section
                    if line.startswith('#') and line.endswith('Coeffs\n'):
                        key = line.lstrip('# ')
                        params[key] = []
                    # strip off things like pair_coeff, bond_coeff, etc
                    elif '_coeff' in first_word and key:
                        params[key].append(line.lstrip(first_word))

            # store data from .data file
            final_file_before_coeffs = []
            final_file_after_coeffs = []
            before_coeffs = True
            with open(f'{tmp_file_prefilx}.data', 'r') as lines:
                for line in lines:
                    if line == 'Atoms\n':
                        before_coeffs = False
                    if before_coeffs:
                        final_file_before_coeffs.append(line)
                    else:
                        final_file_after_coeffs.append(line)

            # there are double empty lines at the end of EMC data file
            final_file_after_coeffs = final_file_after_coeffs[:-1]

            # combine data and parameters into the final data file
            with open(data_fname, 'w') as f:
                for line in final_file_before_coeffs:
                    f.write(line)
                for param, param_lines in params.items():
                    f.write(param)
                    f.write('\n')
                    for line in param_lines:
                        if param == 'Pair Coeffs\n':
                            # remove the extra type id in the line
                            # ex: 1 1    0.05    4.00 -> 1    0.05    4.00
                            stripped = line.split(maxsplit=1)[1]
                            first, second = stripped.split(maxsplit=1)
                            line = (f'{first:>8}    {second}')
                        else:
                            for coeff in EMC_COEFF_EXCLUSIONS:
                                line = line.replace(f' {coeff} ', ' ')
                        f.write(line)
                    f.write('\n')
                for line in final_file_after_coeffs:
                    f.write(line)

            Pmdlogging.info(f'System file - {data_fname} '
                            f'successfully created in {output_dir}')

        finally:
            # Clean up all EMC generated files except for the data file
            if cleanup:
                fnames = ['build.emc']
                for ext in EMC_EXTS:
                    fnames.append(f'{tmp_file_prefilx}.{ext}')
                for fname in fnames:
                    if os.path.isfile(fname):
                        os.remove(fname)

            os.chdir(previous_dir)

    @build_dir
    def write_data(self, output_dir: str, smiles: str, density: float,
                   natoms_total: int, length: int, nchains: int,
                   end_cap_smiles: str, data_fname: str,
                   cleanup: bool) -> None:

        smiles = self._remove_brackets_around_asterisks(smiles)
        end_cap_smiles = self._remove_brackets_around_asterisks(end_cap_smiles)

        tmp_file_prefilx = 'system'
        # Write .esh file required to run EMC
        esh_file = os.path.join(output_dir, f'{tmp_file_prefilx}.esh')
        with open(esh_file, 'w') as f:
            f.write('# EMC input file generated by PMD\n')
            f.write('ITEM OPTIONS\n')
            f.write('replace true\n')
            f.write(f'field {self._force_field}\n')
            f.write(f'density {density}\n')
            f.write(f'ntotal {natoms_total}\n')
            f.write('ITEM END\n')
            f.write('\n')
            f.write('ITEM GROUPS\n')
            f.write(f'RU {smiles},1,RU:2\n')
            f.write(f'terminator {end_cap_smiles},1,RU:1,1,RU:2\n')
            f.write('ITEM END\n')
            f.write('\n')
            f.write('ITEM CLUSTERS\n')
            f.write('poly alternate 1\n')
            f.write('ITEM END\n')
            f.write('\n')
            f.write('ITEM POLYMERS\n')
            f.write('poly\n')
            f.write(f'1 RU,{length},terminator,2\n')
            f.write('ITEM END\n')

        self._run_emc(tmp_file_prefilx, output_dir, data_fname, cleanup)

    @build_dir
    def write_solvent_data(self, output_dir: str, smiles: str,
                           solvent_smiles: str, density: float,
                           natoms_total: int, length: int, nsolvents: int,
                           nchains: int, end_cap_smiles: str, data_fname: str,
                           cleanup: bool) -> None:

        smiles = self._remove_brackets_around_asterisks(smiles)
        end_cap_smiles = self._remove_brackets_around_asterisks(end_cap_smiles)

        tmp_file_prefilx = 'solventsystem'
        # Write .esh file required to run EMC
        esh_file = os.path.join(output_dir, f'{tmp_file_prefilx}.esh')
        with open(esh_file, 'w') as f:
            f.write('# EMC input file generated by PMD\n')
            f.write('ITEM OPTIONS\n')
            f.write('replace true\n')
            f.write(f'field {self._force_field}\n')
            f.write(f'density {density}\n')
            f.write(f'ntotal {natoms_total}\n')
            f.write('ITEM END\n')
            f.write('\n')
            f.write('ITEM GROUPS\n')
            f.write(f'solvent {solvent_smiles}\n')
            f.write(f'RU {smiles},1,RU:2\n')
            f.write(f'terminator {end_cap_smiles},1,RU:1,1,RU:2\n')
            f.write('ITEM END\n')
            f.write('\n')
            f.write('ITEM CLUSTERS\n')
            f.write(f'solvent solvent {nsolvents/(nsolvents + nchains)}\n')
            f.write(f'poly alternate {nchains/(nsolvents + nchains)}\n')
            f.write('ITEM END\n')
            f.write('\n')
            f.write('ITEM POLYMERS\n')
            f.write('poly\n')
            f.write(f'1 RU,{length},terminator,2\n')
            f.write('ITEM END\n')

        self._run_emc(tmp_file_prefilx, output_dir, data_fname, cleanup)

    def write_functional_form(self, f: TextIOWrapper) -> None:
        if self._force_field.startswith('opls'):
            f.write(f'{"pair_style":<15} lj/cut/coul/long 9.5 9.5\n')
            f.write(f'{"pair_modify":<15} mix geometric tail yes\n')
            f.write(f'{"kspace_style":<15} pppm/cg 1e-4\n')
            f.write(f'{"bond_style":<15} harmonic\n')
            f.write(f'{"angle_style":<15} harmonic\n')
            f.write(f'{"dihedral_style":<15} multi/harmonic\n')
            f.write(f'{"improper_style":<15} harmonic\n')
            f.write(f'{"special_bonds":<15} lj/coul 0.0 0.0 0.5\n')

        elif self._force_field == 'pcff':
            f.write(f'{"pair_style":<15} lj/class2/coul/long 9.5 9.5\n')
            f.write(f'{"pair_modify":<15} mix sixthpower tail yes\n')
            f.write(f'{"kspace_style":<15} pppm/cg 1e-4\n')
            f.write(f'{"bond_style":<15} class2\n')
            f.write(f'{"angle_style":<15} class2\n')
            f.write(f'{"dihedral_style":<15} class2\n')
            f.write(f'{"improper_style":<15} class2\n')
            f.write(f'{"special_bonds":<15} lj/coul 0 0 1\n')

        elif self._force_field == 'trappe':
            f.write(f'{"pair_style":<15} lj/cut/coul/long 14.0\n')
            f.write(f'{"pair_modify":<15} mix arithmetic tail yes\n')
            f.write(f'{"kspace_style":<15} pppm/cg 1e-4\n')
            f.write(f'{"bond_style":<15} harmonic\n')
            f.write(f'{"angle_style":<15} harmonic\n')
            f.write(f'{"dihedral_style":<15} multi/harmonic\n')
            f.write(f'{"improper_style":<15} harmonic\n')
            f.write(f'{"special_bonds":<15} lj 0 0 0 coul 0 0 0.5\n')


class PSP(Builder):
    '''Object to perform system structure generation using
    [PSP](https://github.com/Ramprasad-Group/PSP): Polymer Structure Predictor
    package. This object should be used as input argument of `System` or
    `Lammps` objects

    Attributes:
        force_field (str): Force field, options are `"opls-lbcc"`,
            `"opls-cm1a"`, `"gaff2-gasteiger"`, and `"gaff2-am1bcc"`

        packmol_nloop (int): Maximum number of optimization loops of Packmol
            (PSP uses Packmol to pack molecules into a box); default: None

        packmol_precision (float): Packmol avoids atom overlaps by ensuring
            a 2.0 Angs atom distance, this parameter determines how close the
            solution must be to the desired distances to be considered correct
            ; default: None
    '''

    def __init__(self,
                 force_field: str,
                 packmol_nloop: Optional[int] = None,
                 packmol_precision: Optional[float] = None) -> None:

        super().__init__(force_field, PSP_FORCE_FIELD_OPTIONS)
        self._packmol_nloop = packmol_nloop
        self._packmol_precision = packmol_precision

    @staticmethod
    def _add_brackets_to_asterisks(smiles: str) -> str:
        stars_no_bracket = re.findall(r'(?<!\[)\*(?!\])', smiles)
        if len(stars_no_bracket) == 2:
            smiles = smiles.replace("*", "[*]")
        return smiles

    def _is_opls_force_field(self) -> bool:
        return self._force_field.startswith('opls')

    def _run_psp(self, input_data: dict, density: float, data_fname: str,
                 output_dir: str, cleanup: bool) -> None:
        try:
            import psp.AmorphousBuilder as ab
        except ImportError:
            raise ImportError('Please install PSP to use PSP builder')

        Pmdlogging.info(
            'PSP builder: Creating the system, this may take a while...')
        try:
            with HiddenPrints():
                amor = ab.Builder(pd.DataFrame(data=input_data),
                                  density=density,
                                  outdir=output_dir)
                # adjust packmol parameters if specified
                if self._packmol_nloop or self._packmol_precision:
                    amor.set_packmol_params = ab.packmol_params(
                        nloop=self._packmol_nloop,
                        precision=self._packmol_precision)
                amor.Build()

                if self._is_opls_force_field():
                    amor.get_opls(
                        output_fname=data_fname,
                        lbcc_charges=self._force_field.endswith('lbcc'))
                else:
                    amor.get_gaff2(
                        output_fname=data_fname,
                        atom_typing='antechamber',
                        am1bcc_charges=self._force_field.endswith('am1bcc'),
                        swap_dict={
                            'ns': 'n',
                            'nt': 'n',
                            'nv': 'nh'
                        })
            Pmdlogging.info(f'System file - {data_fname} '
                            f'successfully created in {output_dir}')
        finally:
            if cleanup:
                force_field_dname = 'ligpargen' if self._is_opls_force_field(
                ) else 'pysimm'
                dnames = ['chain_models', 'packmol']
                dnames.append(force_field_dname)
                for dir in dnames:
                    dir_path = os.path.join(output_dir, dir)
                    if os.path.isdir(dir_path):
                        shutil.rmtree(dir_path)

                fnames = ['amor_model.data', 'amor_model.vasp']
                for file in fnames:
                    file_path = os.path.join(output_dir, file)
                    if os.path.isfile(file_path):
                        os.remove(file_path)

                fnames = ['output_MB.csv', 'molecules.csv']
                for file in fnames:
                    if os.path.isfile(file):
                        os.remove(file)

    def write_data(self, output_dir: str, smiles: str, density: float,
                   natoms_total: int, length: int, nchains: int,
                   end_cap_smiles: str, data_fname: str,
                   cleanup: bool) -> None:

        smiles = self._add_brackets_to_asterisks(smiles)
        end_cap_smiles = self._add_brackets_to_asterisks(end_cap_smiles)
        input_data = {
            'ID': ['Poly'],
            'smiles': [smiles],
            'Tunits': [length],
            'Num': [nchains],
            'Loop': [False],
            'LeftCap': [end_cap_smiles],
            'RightCap': [end_cap_smiles]
        }
        self._run_psp(input_data, density, data_fname, output_dir, cleanup)

    def write_solvent_data(self, output_dir: str, smiles: str,
                           solvent_smiles: str, density: float,
                           natoms_total: int, length: int, nsolvents: int,
                           nchains: int, end_cap_smiles: str, data_fname: str,
                           cleanup: bool) -> None:

        smiles = self._add_brackets_to_asterisks(smiles)
        end_cap_smiles = self._add_brackets_to_asterisks(end_cap_smiles)
        input_data = {
            'ID': ['Sol', 'Poly'],
            'smiles': [solvent_smiles, smiles],
            'Tunits': [1, length],
            'Num': [nsolvents, nchains],
            'Loop': [False, False],
            'LeftCap': [np.nan, end_cap_smiles],
            'RightCap': [np.nan, end_cap_smiles]
        }
        self._run_psp(input_data, density, data_fname, output_dir, cleanup)

    def write_functional_form(self, f: TextIOWrapper) -> None:

        if self._force_field.startswith('opls'):
            f.write(f'{"pair_style":<15} lj/cut/coul/long 9.0\n')
            f.write(f'{"pair_modify":<15} mix geometric tail yes\n')
            f.write(f'{"kspace_style":<15} pppm 1e-4\n')
            f.write(f'{"bond_style":<15} harmonic\n')
            f.write(f'{"angle_style":<15} harmonic\n')
            f.write(f'{"dihedral_style":<15} opls\n')
            f.write(f'{"improper_style":<15} cvff\n')
            f.write(f'{"special_bonds":<15} lj/coul 0.0 0.0 0.5\n')

        elif self._force_field.startswith('gaff2'):
            f.write(f'{"pair_style":<15} lj/cut/coul/long 12.0 12.0\n')
            f.write(f'{"pair_modify":<15} mix arithmetic\n')
            f.write(f'{"kspace_style":<15} pppm 1e-4\n')
            f.write(f'{"bond_style":<15} harmonic\n')
            f.write(f'{"angle_style":<15} harmonic\n')
            f.write(f'{"dihedral_style":<15} fourier\n')
            f.write(f'{"improper_style":<15} cvff\n')
            f.write(f'{"special_bonds":<15} amber\n')
