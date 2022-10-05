from io import TextIOWrapper
from typing import Optional

from pmd.util.Log import Pmdlogging


class Procedure():

    def __init__(self,
                 duration: int,
                 dump_fname: str,
                 dump_every: int = 10000,
                 dump_image: bool = False,
                 reset_timestep_before_run: bool = False):
        self._duration = duration
        self._dump_fname = dump_fname
        self._dump_every = dump_every
        self._dump_image = dump_image
        self._reset_timestep_before_run = reset_timestep_before_run

    def __repr__(self) -> str:
        return type(self).__name__

    def write_lammps(self, f: TextIOWrapper):
        raise NotImplementedError

    def write_before_run(self, f: TextIOWrapper):
        f.write(f'### {self}\n')
        if self._reset_timestep_before_run:
            f.write(f'{"reset_timestep":<15} 0\n')
            f.write('\n')
        f.write(f'{"dump":<15} dump_{self} all custom {self._dump_every} '
                f'{self._dump_fname} id mol type q xs ys zs ix iy iz\n')
        if self._dump_image:
            f.write(f'{"dump":<15} dump_image all image {self._duration} '
                    f'{self}.*.jpg type type\n')
        f.write(f'{"restart":<15} {self._duration} {self}.restart\n')
        f.write('\n')

    def write_after_run(self, f: TextIOWrapper):
        f.write(f'{"undump":<15} dump_{self}\n')
        if self._dump_image:
            f.write(f'{"undump":<15} dump_image\n')


class Minimization(Procedure):
    '''Perform an energy minimization of the system, by iteratively adjusting
    atom coordinates. Iterations are terminated when one of the stopping
    criteria is satisfied. At that point the configuration will hopefully be in
    local potential energy minimum.

    Attributes:
        min_style (str): Minimization algorithm, see
            [here](https://docs.lammps.org/min_style.html) for all options
            ; default: `"cg"`

        etol (float): Stopping tolerance for energy (unitless); default:
           `10**(-8)`

        ftol (float): Stopping tolerance for force (force units); default:
            `10**(-10)`

        maxiter (int): Max iterations of minimizer; default: `10**7`

        maxeval (int): Max number of force/energy evaluations; default: `10**9`
    '''

    def __init__(self,
                 min_style: str = 'cg',
                 etol: float = 10**(-8),
                 ftol: float = 10**(-10),
                 maxiter: int = 10**7,
                 maxeval: int = 10**9):
        self._min_style = min_style
        self._etol = etol
        self._ftol = ftol
        self._maxiter = maxiter
        self._maxeval = maxeval

    def write_lammps(self, f: TextIOWrapper):
        f.write('### Minimization\n')
        f.write(f'{"min_style":<15} {self._min_style}\n')
        f.write(f'{"minimize":<15} {self._etol} {self._ftol} '
                f'{self._maxiter} {self._maxeval}\n')

    def write_before_run(self, f: TextIOWrapper):
        pass

    def write_after_run(self, f: TextIOWrapper):
        pass


class Equilibration(Procedure):
    '''Perform a 21-step amorphous polymer equilibration process. Ref: Abbott,
    Hart, and Colina, Theoretical Chemistry Accounts, 132(3), 1-19, 2013.

    Attributes:
        Teq (float): Target equilibration temperature; default: `300`

        Peq (float): Target equilibration pressure; default: `1`

        Tmax (float): Maximum temperature during the equilibration; default:
           `600`

        Pmax (float): Maximum pressure during the equilibration; default:
            `50000`

        Tdamp (str): Damping parameter for the thermostat; default:
            `"$(100.0*dt)"`

        Pdamp (str): Damping parameter for the barostat; default:
            `"$(100.0*dt)"`

        nve_limit_start (bool): Whether to start the simulation with a fix
            nve/limit for 10000 timesteps. This avoids simulation failure due
            to a bad initial configuration, see more at
            [here](https://docs.lammps.org/fix_nve_limit.html); default:
            `True`

        dump_fname (str): Name of the dump file; default: `"equil.lammpstrj"`

        dump_every (int): Dump every this many timesteps; default: `10000`

        dump_image (bool): Whether to dump a image file at the end of the run
            ; default: `False`

        reset_timestep_before_run (bool): Whether to reset timestep after the
            procedure; default: `True`
    '''

    def __init__(self,
                 Teq: float = 300,
                 Peq: float = 1,
                 Tmax: float = 600,
                 Pmax: float = 50000,
                 Tdamp: str = '$(100.0*dt)',
                 Pdamp: str = '$(100.0*dt)',
                 nve_limit_start: bool = True,
                 dump_fname: str = 'equil.lammpstrj',
                 dump_every: int = 10000,
                 dump_image: bool = False,
                 reset_timestep_before_run: bool = True):
        self._Teq = Teq
        self._Peq = Peq
        self._Tmax = Tmax
        self._Pmax = Pmax
        self._Tdamp = Tdamp
        self._Pdamp = Pdamp
        self._nve_limit_start = nve_limit_start

        duration = 0
        for i in self._eq_steps:
            duration += i[1]

        super().__init__(duration, dump_fname, dump_every, dump_image,
                         reset_timestep_before_run)

    @property
    def _eq_steps(self):
        return [
            ['nvt', 50000, self._Tmax],
            ['nvt', 50000, self._Teq],
            ['npt', 50000, self._Teq, 0.02 * self._Pmax],
            ['nvt', 50000, self._Tmax],
            ['nvt', 100000, self._Teq],
            ['npt', 50000, self._Teq, 0.6 * self._Pmax],
            ['nvt', 50000, self._Tmax],
            ['nvt', 100000, self._Teq],
            ['npt', 50000, self._Teq, self._Pmax],
            ['nvt', 50000, self._Tmax],
            ['nvt', 100000, self._Teq],
            ['npt', 5000, self._Teq, 0.5 * self._Pmax],
            ['nvt', 5000, self._Tmax],
            ['nvt', 10000, self._Teq],
            ['npt', 5000, self._Teq, 0.1 * self._Pmax],
            ['nvt', 5000, self._Tmax],
            ['nvt', 10000, self._Teq],
            ['npt', 5000, self._Teq, 0.01 * self._Pmax],
            ['nvt', 5000, self._Tmax],
            ['nvt', 10000, self._Teq],
            ['npt', 800000, self._Teq, self._Peq],
        ]

    def write_lammps(self, f: TextIOWrapper):
        if self._nve_limit_start:
            f.write(f'{"fix":<15} fLANGEVIN all langevin '
                    f'{self._Tmax} {self._Tmax} {self._Tdamp} 723853\n')
            f.write(f'{"fix":<15} fNVELIMIT all nve/limit 0.1\n')
            f.write(
                f'{"fix":<15} fMOM all momentum 100 linear 1 1 1 angular\n')
            f.write(f'{"run":<15} 10000\n')
            f.write('\n')
            f.write(f'{"unfix":<15} fLANGEVIN\n')
            f.write(f'{"unfix":<15} fNVELIMIT\n')
            f.write(f'{"unfix":<15} fMOM\n')
            f.write('\n')
            f.write(f'{"reset_timestep":<15} 0\n')
            f.write('\n')

        for n, i in enumerate(self._eq_steps):
            if i[0] == 'nvt':
                f.write(f'{"fix":<15} step{n + 1} all nvt temp '
                        f'{i[2]} {i[2]} {self._Tdamp}\n')
            elif i[0] == 'npt':
                f.write(f'{"fix":<15} step{n + 1} all npt temp {i[2]} {i[2]} '
                        f'{self._Tdamp} iso {i[3]} {i[3]} {self._Pdamp}\n')
            f.write(f'{"run":<15} {i[1]}\n')
            f.write(f'{"unfix":<15} step{n + 1}\n')
            f.write('\n')


class NPT(Procedure):
    '''Perform the simulation under NPT ensemble (via Nose-Hoover thermostat
    and barostat).

    Attributes:
        duration (int): Duration of this NPT procedure (timestep unit)

        Tinit (float): Initial temperature

        Tfinal (float): Final temperature

        Pinit (float): Initial pressure

        Pfinal (float): Final pressure

        Tdamp (str): Damping parameter for the thermostat; default:
            `"$(100.0*dt)"`

        Pdamp (str): Damping parameter for the barostat; default:
            `"$(100.0*dt)"`

        dump_fname (str): Name of the dump file; default: `"npt.lammpstrj"`

        dump_every (int): Dump every this many timesteps; default: `10000`

        dump_image (bool): Whether to dump a image file at the end of the run
            ; default: `False`

        reset_timestep_before_run (bool): Whether to reset timestep after the
            procedure; default: `False`
    '''

    def __init__(self,
                 duration: int,
                 Tinit: float,
                 Tfinal: float,
                 Pinit: float,
                 Pfinal: float,
                 Tdamp: str = '$(100.0*dt)',
                 Pdamp: str = '$(100.0*dt)',
                 dump_fname: str = 'npt.lammpstrj',
                 dump_every: int = 10000,
                 dump_image: bool = False,
                 reset_timestep_before_run: bool = False):
        self._Tinit = Tinit
        self._Tfinal = Tfinal
        self._Pinit = Pinit
        self._Pfinal = Pfinal
        self._Tdamp = Tdamp
        self._Pdamp = Pdamp

        super().__init__(duration, dump_fname, dump_every, dump_image,
                         reset_timestep_before_run)

    def write_lammps(self, f: TextIOWrapper):
        f.write(
            f'{"fix":<15} fNPT all npt temp {self._Tinit} {self._Tfinal} '
            f'{self._Tdamp} iso {self._Pinit} {self._Pfinal} {self._Pdamp}\n')
        f.write(f'{"run":<15} {self._duration}\n')
        f.write(f'{"unfix":<15} fNPT\n')
        f.write('\n')


class NVT(Procedure):
    '''Perform the simulation under NVT ensemble (via Nose-Hoover thermostat).

    Attributes:
        duration (int): Duration of this NVT procedure (timestep unit)

        Tinit (float): Initial temperature

        Tfinal (float): Final temperature

        Tdamp (str): Damping parameter for thermostats; default:
            `"$(100.0*dt)"`

        dump_fname (str): Name of the dump file; default: `"nvt.lammpstrj"`

        dump_every (int): Dump every this many timesteps; default: `10000`

        dump_image (bool): Whether to dump a image file at the end of the run
            ; default: `False`

        reset_timestep_before_run (bool): Whether to reset timestep after the
            procedure; default: `False`
    '''

    def __init__(self,
                 duration: int,
                 Tinit: float,
                 Tfinal: float,
                 Tdamp: str = '$(100.0*dt)',
                 dump_fname: str = 'nvt.lammpstrj',
                 dump_every: int = 10000,
                 dump_image: bool = False,
                 reset_timestep_before_run: bool = False):
        self._Tinit = Tinit
        self._Tfinal = Tfinal
        self._Tdamp = Tdamp

        super().__init__(duration, dump_fname, dump_every, dump_image,
                         reset_timestep_before_run)

    def write_lammps(self, f: TextIOWrapper):
        f.write(f'{"fix":<15} fNVT all nvt temp {self._Tinit} '
                f'{self._Tfinal} {self._Tdamp}\n')
        f.write(f'{"run":<15} {self._duration}\n')
        f.write(f'{"unfix":<15} fNVT\n')
        f.write('\n')


class MSDMeasurement(Procedure):
    '''Perform mean-squared displacement measurement for the specified group
    of atmos/molecules.

    Attributes:
        duration (int): Duration of the NVT ensemble for MSD measurement
            (timestep unit)

        T (float): Temperature

        group (str): The group of atoms that will be considered for MSD
            calculation. This has to be a string that matches the syntax of
            [group](https://docs.lammps.org/group.html) LAMMPS command
            (e.g. `"molecule <=50"`, `"type 1 2"`, etc

        create_block_every (int): The time interval that new MSD calculation
            starting point will be created (e.g. for a 1000 fs run, a
            `create_block_every` value of 100fs would result in 10 blocks with
            10 different MSD starting point and length) ; default: `None`

        result_folder_name (str): The name of the folder that PMD creates and
            put result files in; default: `"result"`

        Tdamp (str): Damping parameter for thermostats; default:
            `"$(100.0*dt)"`

        dump_fname (str): Name of the dump file; default: `"nvt.lammpstrj"`

        dump_every (int): Dump every this many timesteps; default: `10000`

        dump_image (bool): Whether to dump a image file at the end of the run
            ; default: `False`

        reset_timestep_before_run (bool): Whether to reset timestep after the
            procedure; default: `False`
    '''

    def __init__(self,
                 duration: int,
                 T: float,
                 group: str,
                 create_block_every: Optional[int] = None,
                 result_folder_name: str = 'result',
                 Tdamp: str = '$(100.0*dt)',
                 dump_fname: str = 'MSD_measurement.lammpstrj',
                 dump_every: int = 10000,
                 dump_image: bool = False,
                 reset_timestep_before_run: bool = False):

        if duration % create_block_every != 0:
            raise ValueError('The duration has to be divisible by the '
                             'create_block_every')

        self._T = T
        self._group = group
        self._create_block_every = create_block_every
        self._Tdamp = Tdamp
        self._result_folder_name = result_folder_name

        super().__init__(duration, dump_fname, dump_every, dump_image,
                         reset_timestep_before_run)

    def write_lammps(self, f: TextIOWrapper):
        f.write(f'{"fix":<15} fNVT all nvt '
                f'temp {self._T} {self._T} {self._Tdamp}\n')
        f.write('\n')

        msd_group_id = 'msdgroup'
        mol_chunk_id = 'molchunk'
        msd_chunk_id = 'msdchunk'
        f.write(f'{"shell":<15} mkdir {self._result_folder_name}\n')
        f.write(f'{"group":<15} {msd_group_id} {self._group}\n')
        f.write(f'{"compute":<15} {mol_chunk_id} {msd_group_id} '
                f'chunk/atom molecule\n')
        f.write('\n')

        if self._create_block_every:
            nblock = int(self._duration / self._create_block_every)
        else:
            nblock = 1
        for block in range(nblock):
            start = block * self._create_block_every
            f.write(f'##### MSDMeasurement block {block}\n')
            f.write(f'{"compute":<15} {msd_chunk_id}{block} {msd_group_id} '
                    f'msd/chunk {mol_chunk_id}\n')
            f.write(f'{"variable":<15} ave{msd_chunk_id}{block} equal '
                    f'ave(c_{msd_chunk_id}{block}[4])\n')
            f.write(
                f'{"fix":<15} fMSD{block} {msd_group_id} ave/time '
                f'1 1 10000 v_ave{msd_chunk_id}{block} start {start} file '
                f'{self._result_folder_name}/msd_{start}_{self._duration}.txt'
                f'\n')
            f.write(f'{"run":<15} {self._create_block_every}\n')
            f.write('\n')

        f.write('\n')
        f.write(f'{"unfix":<15} fNVT\n')
        for block in range(nblock):
            f.write(f'{"unfix":<15} fMSD{block}\n')
        f.write('\n')


class TgMeasurement(Procedure):
    '''Perform glass transition temperature measurement of the system,
    by iteratively cooling the system and equilibrate.

    Attributes:
        Tinit (float): Initial temperature of the cooling process; default:
            `500`

        Tfinal (float): Final temperature of the cooling process; default:
            `100`

        Tinterval (float): Temperature interval of the cooling process
            ; default: `20`

        step_duration (int): Duration of each temperature step (timestep unit)
            ; default: `1000000`

        pressure (float): Pressure during the cooling process; default: `1`

        Tdamp (str): Damping parameter for the thermostat; default:
            `"$(100.0*dt)"`

        Pdamp (str): Damping parameter for the barostat; default:
            `"$(100.0*dt)"`

        dump_fname (str): Name of the dump file; default:
            `"Tg_measurement.lammpstrj"`

        dump_every (int): Dump every this many timesteps; default: `10000`

        dump_image (bool): Whether to dump a image file at the end of the run
            ; default: `False`

        result_fname (str): Name of the result file; default:
            `"temp_vs_density.txt"`
    '''

    def __init__(self,
                 Tinit: float = 500,
                 Tfinal: float = 100,
                 Tinterval: float = 20,
                 step_duration: int = 1000000,
                 pressure: float = 1,
                 Tdamp: str = '$(100.0*dt)',
                 Pdamp: str = '$(100.0*dt)',
                 result_fname: str = 'temp_vs_density.txt',
                 dump_fname: str = 'Tg_measurement.lammpstrj',
                 dump_every: int = 10000,
                 dump_image: bool = False,
                 reset_timestep_before_run: bool = False):
        if Tinit < Tfinal:
            Pmdlogging.warning(
                f'Tg measurement usually is done through a cooling process,'
                f' but your Tinit: {Tinit} is lower than Tfinal: {Tfinal}')
        if (Tinit - Tfinal) % Tinterval != 0:
            Pmdlogging.warning(
                'Your Tinterval is not a factor of Tinit-Tfinal')
        self._Tinit = Tinit
        self._Tfinal = Tfinal
        self._Tinterval = Tinterval
        self._step_duration = step_duration
        self._pressure = pressure
        self._Tdamp = Tdamp
        self._Pdamp = Pdamp
        self._result_fname = result_fname

        # total duration of the procedure
        duration = self._nsteps * step_duration

        super().__init__(duration, dump_fname, dump_every, dump_image,
                         reset_timestep_before_run)

    @property
    def _nsteps(self) -> int:
        return int((self._Tinit - self._Tfinal) / self._Tinterval + 1)

    def write_lammps(self, f: TextIOWrapper):
        f.write(f'{"variable":<15} Rho equal density\n')
        f.write(f'{"variable":<15} Temp equal temp\n')
        f.write(
            f'{"fix":<15} fDENS all ave/time '
            f'{int(self._step_duration / 100 / 4)} {100} '
            f'{self._step_duration} v_Temp v_Rho file {self._result_fname}\n')
        f.write('\n')

        f.write(f'{"label":<15} loop\n')
        f.write(f'{"variable":<15} a loop {self._nsteps}\n')
        f.write(f'{"variable":<15} b equal '
                f'{self._Tinit}-{self._Tinterval}*($a-1)\n')
        f.write(f'{"fix":<15} fNPT all npt temp $b $b {self._Tdamp} iso '
                f'{self._pressure} {self._pressure} {self._Pdamp}\n')
        f.write(f'{"run":<15} {self._step_duration}\n')
        f.write(f'{"unfix":<15} fNPT\n')
        f.write(f'{"next":<15} a\n')
        f.write(f'{"jump":<15} SELF loop\n')
        f.write(f'{"variable":<15} a delete\n')
        f.write('\n')


class TensileDeformation(Procedure):
    '''Perform a uniaxial tensile deformation in the x direction.
    This can be used to calculate modulus and tensile strengths.

    Attributes:
        duration (int): Duration of the deformation procedure (timestep unit)

        erate (float): Engineering strain rate. The units of the specified
            strain rate are 1/time

        T (float): Temperature

        P (float): Pressure

        Tdamp (str): Damping parameter for thermostats; default:
            `"$(100.0*dt)"`

        Pdamp (str): Damping parameter for thermostats; default:
            `"$(100.0*dt)"`

        print_every (int): Print result to the result file every this many
            timesteps; default: `1000`

        dump_fname (str): Name of the dump file; default:
            `"tensile_deformation.lammpstrj"`

        dump_every (int): Dump every this many timesteps; default: `10000`

        dump_image (bool): Whether to dump a image file at the end of the run
            ; default: `False`

        reset_timestep_before_run (bool): Whether to reset timestep after the
             procedure; default: `False`

        result_fname (str): Name of the result file; default:
            `"stress_vs_strain.txt"`
    '''

    def __init__(self,
                 duration: int,
                 erate: float,
                 T: float,
                 P: float,
                 Tdamp: str = '$(100.0*dt)',
                 Pdamp: str = '$(100.0*dt)',
                 print_every: int = 1000,
                 result_fname: str = 'stress_vs_strain.txt',
                 dump_fname: str = 'tensile_deformation.lammpstrj',
                 dump_every: int = 10000,
                 dump_image: bool = False,
                 reset_timestep_before_run: bool = False):
        self._erate = erate
        self._T = T
        self._P = P
        self._Tdamp = Tdamp
        self._Pdamp = Pdamp
        self._print_every = print_every
        self._result_fname = result_fname

        super().__init__(duration, dump_fname, dump_every, dump_image,
                         reset_timestep_before_run)

    def write_lammps(self, f: TextIOWrapper):
        f.write(f'{"fix":<15} fNPT all npt '
                f'temp {self._T} {self._T} {self._Tdamp} '
                f'y {self._P} {self._P} {self._Pdamp} '
                f'z {self._P} {self._P} {self._Pdamp}\n')
        f.write(f'{"fix":<15} fDEFORM all deform 1 x erate {self._erate} '
                f'units box remap x\n')

        f.write('\n')
        f.write(f'{"variable":<15} tmp equal "lx"\n')
        f.write(f'{"variable":<15} L0 equal {"${tmp}"}\n')
        f.write(f'{"variable":<15} strain equal "(lx - v_L0)/v_L0"\n')
        f.write(f'{"variable":<15} p1 equal "v_strain"\n')

        # TODO: this assume the LAMMPS units is real, make it more flexible
        # Output strain and stress info to file for units real,
        # pressure is in [atm] = 0.101325 [MPa]
        # and p2, p3, p4 are in MPa
        f.write(f'{"variable":<15} p2 equal "-pxx*0.101325" '
                '# convert stress unit from atm to MPa\n')
        f.write(f'{"variable":<15} p3 equal "-pyy*0.101325"\n')
        f.write(f'{"variable":<15} p4 equal "-pzz*0.101325"\n')
        f.write(f'{"fix":<15} fAVETIME all ave/time 100 {self._print_every} '
                f'{100*self._print_every} v_p1 v_p2 ave one file '
                f'{self._result_fname}\n')
        # override the default thermo_style
        f.write(f'{"thermo_style":<15} custom step temp density vol v_strain '
                'v_p2 v_p3 v_p4 press ke pe\n')
        f.write(f'{"run":<15} {self._duration}\n')
        f.write('\n')
        f.write(f'{"unfix":<15} fNPT\n')
        f.write(f'{"unfix":<15} fDEFORM\n')
        f.write(f'{"unfix":<15} fAVETIME\n')
        f.write('\n')


class ShearDeformation(Procedure):
    '''Perform a shear deformation in the x-y plane. This can be used
    to calculate shear viscosity.

    Attributes:
        duration (int): Duration of the deformation procedure (timestep unit)

        shear_rate (float): Shear rate [1/s] (engineering strain rate
            in LAMMPS, see [here](https://docs.lammps.org/fix_deform.html))

        T (float): Temperature

        Tdamp (str): Damping parameter for thermostats; default:
            `"$(100.0*dt)"`

        calculate_every (int): Calculate result every this many
            timesteps; default: `100000`

        dump_fname (str): Name of the dump file; default:
            `"shear_deformation.lammpstrj"`

        dump_every (int): Dump every this many timesteps; default: `10000`

        dump_image (bool): Whether to dump a image file at the end of the run
            ; default: `False`

        reset_timestep_before_run (bool): Whether to reset timestep after the
            procedure; default: `False`

        result_fname (str): Name of the result file, viscosity will be dumped
            out to this file in the unit of [Pa s]; default: `"viscosity.txt"`
    '''

    def __init__(self,
                 duration: int,
                 shear_rate: float,
                 T: float,
                 Tdamp: str = '$(100.0*dt)',
                 calculate_every: int = 100000,
                 result_fname: str = 'viscosity.txt',
                 dump_fname: str = 'shear_deformation.lammpstrj',
                 dump_every: int = 10000,
                 dump_image: bool = False,
                 reset_timestep_before_run: bool = False):
        self._shear_rate = shear_rate
        self._T = T
        self._Tdamp = Tdamp
        self._calculate_every = calculate_every
        self._result_fname = result_fname

        super().__init__(duration, dump_fname, dump_every, dump_image,
                         reset_timestep_before_run)

    # dump instructions have to go after change_box command
    def write_before_run(self, f: TextIOWrapper):
        pass

    def write_lammps(self, f: TextIOWrapper):
        f.write(f'### {self}\n')
        if self._reset_timestep_before_run:
            f.write(f'{"reset_timestep":<15} 0\n')
        f.write('\n')

        f.write(f'{"change_box":<15} all triclinic\n')
        f.write(f'{"kspace_style":<15} pppm 1e-4  '
                '# must redefine pppm after changing to triclinic\n')
        f.write('\n')

        f.write(f'{"dump":<15} dump_{self} all custom {self._dump_every} '
                f'{self._dump_fname} id mol type q xs ys zs ix iy iz\n')
        if self._dump_image:
            f.write(f'{"dump":<15} dump_image all image {self._duration} '
                    f'{self}.*.jpg type type\n')
        f.write(f'{"restart":<15} {self._duration} {self}.restart\n')
        f.write('\n')

        f.write(f'{"variable":<15} srate_in_s equal {self._shear_rate}  '
                '# shear rate [1/s]\n')
        f.write(f'{"variable":<15} srate equal '
                '${srate_in_s}/1e15  '
                '# convert shear rate unit from 1/s to 1/fs\n')
        f.write('\n')

        f.write(f'{"fix":<15} fNVTSLLOD all nvt/sllod '
                f'temp {self._T} {self._T} {self._Tdamp}\n')
        f.write(f'{"fix":<15} fDEFORM all deform 1 xy erate '
                '${srate} remap v\n')
        f.write('\n')

        temp_deform_id = 'tdeform'
        press_deform_id = 'pdeform'
        f.write(f'{"compute":<15} {temp_deform_id} all temp/deform  '
                '# calculate temperautre by subtracting out a '
                'streaming velocity induced by deformation\n')
        f.write(f'{"compute":<15} {press_deform_id} all '
                f'pressure {temp_deform_id}\n')
        f.write(f'{"thermo_modify":<15} temp {temp_deform_id}\n')
        f.write(f'{"thermo_modify":<15} press {press_deform_id}\n')
        f.write(f'{"variable":<15} stress equal '
                f'(-1)*c_{press_deform_id}[4] # -pxy\n')

        # TODO: this assume the LAMMPS units is real, make it more flexible
        # Caclulate shear viscosity [Pa s]
        f.write(f'{"variable":<15} visc equal 1.01325e-10*v_stress/(v_srate) '
                ' # shear viscosity [Pa s]; first term is unit converter\n')
        f.write(f'{"fix":<15} fVISC all '
                f'ave/time 100 1000 {self._calculate_every} v_visc '
                f'ave one file {self._result_fname}\n')
        f.write(f'{"run":<15} {self._duration}\n')
        f.write('\n')
        f.write(f'{"unfix":<15} fNVTSLLOD\n')
        f.write(f'{"unfix":<15} fDEFORM\n')
        f.write(f'{"unfix":<15} fVISC\n')
        f.write('\n')


class HeatFluxMeasurement(Procedure):
    '''Perform a heat flux measurement to calculate the thermal conductivity
    using the equilibrium Green-Kubo formalism.

    Attributes:
        duration (int): Duration of the deformation procedure (timestep unit)

        T (float): Temperature

        Tdamp (str): Damping parameter for thermostats; default:
            `"$(100.0*dt)"`

        dump_fname (str): Name of the dump file; default:
            `"heatflux_measurement.lammpstrj"`

        dump_every (int): Dump every this many timesteps; default: `10000`

        dump_image (bool): Whether to dump a image file at the end of the run
            ; default: `False`

        reset_timestep_before_run (bool): Whether to reset timestep after the
            procedure; default: `False`

        result_fname (str): Name of the result file; default: `"J0Jt.txt"`
    '''

    def __init__(self,
                 duration: int,
                 T: float,
                 Tdamp: str = '$(100.0*dt)',
                 result_fname: str = 'J0Jt.txt',
                 dump_fname: str = 'heatflux_measurement.lammpstrj',
                 dump_every: int = 10000,
                 dump_image: bool = False,
                 reset_timestep_before_run: bool = False):
        self._T = T
        self._Tdamp = Tdamp
        self._result_fname = result_fname

        super().__init__(duration, dump_fname, dump_every, dump_image,
                         reset_timestep_before_run)

    def write_lammps(self, f: TextIOWrapper):
        f.write(f'{"variable":<15} T equal {self._T}\n')
        f.write(f'{"variable":<15} V equal vol\n')
        f.write(f'{"variable":<15} p equal 200\n')
        f.write(f'{"variable":<15} s equal 10\n')
        f.write(f'{"variable":<15} d equal $p*$s\n')
        f.write(f'{"thermo":<15} $d\n')
        f.write('\n')

        f.write('# convert from LAMMPS real units to SI\n')
        f.write(f'{"variable":<15} kB equal 1.3806504e-23  '
                '# [J/K] Boltzmann\n')
        f.write(f'{"variable":<15} kCal2J equal 4186.0/6.02214e23\n')
        f.write(f'{"variable":<15} A2m equal 1.0e-10\n')
        f.write(f'{"variable":<15} fs2s equal 1.0e-15\n')
        f.write(f'{"variable":<15} convert equal '
                '${kCal2J}*${kCal2J}/${fs2s}/${A2m}\n')
        f.write('\n')

        f.write(f'{"fix":<15} fNVT all nvt temp '
                f'{self._T} {self._T} {self._Tdamp} drag 0.2\n')
        f.write('\n')

        f.write(f'{"compute":<15} myKE all ke/atom\n')
        f.write(f'{"compute":<15} myPE all pe/atom\n')
        f.write(f'{"compute":<15} myStress all stress/atom NULL virial\n')
        f.write(f'{"compute":<15} flux all heat/flux myKE myPE myStress\n')
        f.write(f'{"fix":<15} JJ all ave/correlate $s $p $d '
                'c_flux[1] c_flux[2] c_flux[3] type auto file '
                f'{self._result_fname} ave running\n')
        f.write('\n')

        f.write(f'{"variable":<15} scale equal '
                '${convert}/${kB}/$T/$T/$V*$s\n')
        f.write(f'{"variable":<15} '
                'k11 equal trap(f_JJ[3])*${scale}\n')
        f.write(f'{"variable":<15} '
                'k22 equal trap(f_JJ[4])*${scale}\n')
        f.write(f'{"variable":<15} '
                'k33 equal trap(f_JJ[5])*${scale}\n')
        f.write('\n')
        f.write(f'{"run":<15} {self._duration}\n')
        f.write(f'{"variable":<15} k equal (v_k11+v_k22+v_k33)/3.0\n')
        f.write(f'{"print":<15} \"average conductivity: $k[W/mK]\"\n')
        f.write('\n')
        f.write(f'{"unfix":<15} fNVT\n')
        f.write(f'{"unfix":<15} JJ\n')
        f.write('\n')


class RgMeasurement(Procedure):
    '''Perform radius of gyration measurement for the specified group of
    molecules under a NPT ensemble.

    Attributes:
        duration (int): Duration of this procedure (timestep unit)

        T (float): Temperature

        P (float): Pressure

        group (str): The group of atoms that will be considered for MSD
            calculation. This has to be a string that matches the syntax of
            [group](https://docs.lammps.org/group.html) LAMMPS command
            (e.g. `"molecule <=50"`, `"type 1 2"`, etc

        result_fname (str): Name of the result file; default:
            `"Rg_results.txt"`

        calculate_every (int): Calculate result every this many
            timesteps; default: `100000`

        Tdamp (str): Damping parameter for thermostats; default:
            `"$(100.0*dt)"`

        Pdamp (str): Damping parameter for thermostats; default:
            `"$(100.0*dt)"`

        dump_fname (str): Name of the dump file; default:
            `"Rg_Measurement.lammpstrj"`

        dump_every (int): Dump every this many timesteps; default: `10000`

        dump_image (bool): Whether to dump a image file at the end of the run
            ; default: `False`

        reset_timestep_before_run (bool): Whether to reset timestep after the
            procedure; default: `False`
    '''

    def __init__(self,
                 duration: int,
                 T: float,
                 P: float,
                 group: str,
                 result_fname: str = 'Rg_results.txt',
                 calculate_every: int = 100000,
                 Tdamp: str = '$(100.0*dt)',
                 Pdamp: str = '$(100.0*dt)',
                 dump_fname: str = 'Rg_Measurement.lammpstrj',
                 dump_every: int = 10000,
                 dump_image: bool = False,
                 reset_timestep_before_run: bool = False):

        self._T = T
        self._P = P
        self._group = group
        self._result_fname = result_fname
        self._calculate_every = calculate_every
        self._Tdamp = Tdamp
        self._Pdamp = Pdamp

        super().__init__(duration, dump_fname, dump_every, dump_image,
                         reset_timestep_before_run)

    def write_lammps(self, f: TextIOWrapper):
        f.write(f'{"fix":<15} fNPT all npt temp {self._T} {self._T} '
                f'{self._Tdamp} iso {self._P} {self._P} {self._Pdamp}\n')
        f.write('\n')

        rg_group_id = 'rggroup'
        f.write(f'{"group":<15} {rg_group_id} {self._group}\n')
        # TODO: make the mol chunk work for Rg calculation, currently,
        # this only works for the single polymer case
        # mol_chunk_id = 'molchunk'
        # rg_chunk_id = 'rgchunk'
        # f.write(f'{"compute":<15} {mol_chunk_id} {rg_group_id} '
        #         f'chunk/atom molecule\n')
        # f.write('\n')

        # f.write(f'{"compute":<15} {rg_chunk_id} {rg_group_id} '
        #         f'gyration/chunk {mol_chunk_id}  # Rg of each molecule\n')
        # f.write(f'{"variable":<15} Rg equal ave(c_{rg_chunk_id})  '
        #         '# average Rg of all molecules\n')
        f.write(f'{"compute":<15} Rg {rg_group_id} gyration\n')
        f.write(
            f'{"fix":<15} fAVETIME {rg_group_id} ave/time '
            f'100 {int(self._calculate_every/100)} {self._calculate_every} '
            f'c_Rg ave one file {self._result_fname}\n')
        f.write(f'{"run":<15} {self._duration}\n')
        f.write('\n')

        f.write('\n')
        f.write(f'{"unfix":<15} fNPT\n')
        f.write(f'{"unfix":<15} fAVETIME\n')
        f.write('\n')
