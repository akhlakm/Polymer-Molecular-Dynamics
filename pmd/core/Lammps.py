from pmd.core.Procedure import Procedure
from pmd.core.System import System
from pmd.util import Util
from typing import TypeVar

Lammps = TypeVar("Lammps", bound="Lammps")


class Lammps:
    '''Template object to contain LAMMPS initialization settings

    Attributes:
        read_data_from (System): System object that the data file will be read
                                 from; one of this attribute and 
                                 `read_restart_from` has to be provided but not
                                 both (providing both will result in an error);
                                 default: `None`
        
        read_restart_from (Lammps): Lammps object that the last restart file 
                                    created will be read from; one of this 
                                    attribute and `read_data_from` has to be
                                    provided but not both (providing both will
                                    result in an error); default: `None`

        atom_style (str): LAMMPS 
                          [atom_style](https://docs.lammps.org/atom_style.html)
                          to use during simulation; default: `full`

        units (str): LAMMPS [units](https://docs.lammps.org/units.html) to use
                     during simulation; default: `real`
                     
        timestep (float): LAMMPS 
                          [timestep](https://docs.lammps.org/timestep.html) to
                          use during simulation; default: `1 fs`

        neighbor_skin (float): LAMMPS 
                            [neighbor](https://docs.lammps.org/neighbor.html)
                            skin size to use during the simulation; default: 
                            `2.0 Angstrom`
                            
        neighbor_every (int): LAMMPS 
                            [neighbor](https://docs.lammps.org/neighbor.html) 
                            list checking frequency to use during the
                            simulation; default: `1 fs`

        thermo (int): LAMMPS [thermo](https://docs.lammps.org/thermo.html) 
                      to use during simulation; default: `1000 timestep`
            
        lmp_input_fname (str): Name of the LAMMPS input file; default: `lmp.in`
    '''

    def __init__(self,
                 read_data_from: System = None,
                 read_restart_from: Lammps = None,
                 atom_style: str = 'full',
                 units: str = 'real',
                 timestep: int = 1,
                 neighbor_skin: float = 2.0,
                 neighbor_every: int = 1,
                 thermo: int = 1000,
                 lmp_input_fname: str = 'lmp.in'):

        if not read_data_from and not read_restart_from:
            raise ValueError(
                'One of read_data_from and read_restart_from has to be defined'
            )
        elif read_data_from and read_restart_from:
            raise ValueError(
                'Only one of read_data_from and read_restart_from can be defined'
            )

        self._read_data_from = read_data_from
        self._read_restart_from = read_restart_from
        self._atom_style = atom_style
        self._units = units
        self._timestep = timestep
        self._neighbor_skin = neighbor_skin
        self._neighbor_every = neighbor_every
        self._thermo = thermo
        self._lmp_input_fname = lmp_input_fname
        self._procedures = []

    @property
    def lmp_input_fname(self) -> str:
        return self._lmp_input_fname

    def add_procedure(self, procedure: Procedure) -> Lammps:
        '''Method to add simulation procedure
        Parameters:
            procedure (Procedure): One of `minimization`, `equilibration`, 
            `NPT`, `NVT`, and `Tg_measurement`

        Returns:
            Lammps (Lammps): Lammps instance itself (builder design pattern)
        '''

        self._procedures.append(procedure)
        return self

    def write_lammps(self, output_dir: str = '.') -> None:
        '''Method to make LAMMPS input files
        Parameters:
            output_dir (str): Directory for the generated LAMMPS input file
                              ; default: `.`

        Returns:
            None
        '''

        Util.build_dir(output_dir)

        # Write LAMMPS input file
        with open(output_dir + '/' + self._lmp_input_fname, 'w') as f:

            f.write('# LAMMPS input file generated by PMD package\n')
            f.write('\n')

            f.write('### Initialization\n')
            f.write('{:<15} full\n'.format('atom_style'))
            f.write('{:<15} {}\n'.format('units', self._units))
            if self._read_data_from:
                f.write('{:<15} {}\n'.format('read_data',
                                             self._read_data_from.data_fname))
                f.write('\n')
                # TODO: move to ForceField Object?
                force_field = self._read_data_from.force_field
                if (force_field == 'gaff2'):
                    f.write('{:<15} lj/cut/coul/long 12.0 12.0\n'.format(
                        'pair_style'))
                    f.write('{:<15} mix arithmetic\n'.format('pair_modify'))
                    f.write('{:<15} pppm 1e-4\n'.format('kspace_style'))
                    f.write('{:<15} harmonic\n'.format('bond_style'))
                    f.write('{:<15} harmonic\n'.format('angle_style'))
                    f.write('{:<15} fourier\n'.format('dihedral_style'))
                    f.write('{:<15} cvff\n'.format('improper_style'))
                    f.write('{:<15} amber\n'.format('special_bonds'))
                elif (force_field == 'opls'):
                    f.write(
                        '{:<15} lj/cut/coul/long 9.0\n'.format('pair_style'))
                    f.write('{:<15} mix geometric tail yes\n'.format(
                        'pair_modify'))
                    f.write('{:<15} pppm 1e-4\n'.format('kspace_style'))
                    f.write('{:<15} harmonic\n'.format('bond_style'))
                    f.write('{:<15} harmonic\n'.format('angle_style'))
                    f.write('{:<15} opls\n'.format('dihedral_style'))
                    f.write('{:<15} cvff\n'.format('improper_style'))
                    f.write(
                        '{:<15} lj/coul 0.0 0.0 0.5\n'.format('special_bonds'))

            # TODO: add last_restart_fname
            # elif self._read_restart_from:
            # f.write('{:<15} {}\n'.format(
            #     'read_restart',
            #     self._read_restart_from.last_restart_fname))
            f.write('\n')
            f.write('{:<15} {} bin\n'.format('neighbor', self._neighbor_skin))
            f.write('{:<15} delay 0 every {} check yes\n'.format(
                'neigh_modify', self._neighbor_every))
            f.write('\n')
            f.write(
                '{:<15} custom step temp density vol press ke pe ebond evdwl ecoul elong\n'
                .format('thermo_style'))
            f.write('{:<15} {}\n'.format('thermo', self._thermo))
            f.write('{:<15} {}\n'.format('timestep', self._timestep))
            f.write('\n')
            f.write('\n')

            for procedure in self._procedures:
                procedure.write_lammps(f)
