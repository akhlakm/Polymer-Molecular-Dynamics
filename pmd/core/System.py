from typing import Optional, Tuple

from rdkit import Chem
from rdkit.Chem import Descriptors

from pmd.core.Builder import Builder
from pmd.util import Pmdlogging, validate_options

SYSTEM_SIZE_OPTIONS = ('natoms_total', 'nchains_total')
CHAIN_LENGTH_OPTIONS = ('natoms_per_chain', 'mw_per_chain', 'ru_per_chain')


class System:
    '''Template object to contain System initialization settings

    Attributes:
        smiles (str): SMILES string of the polymer (use * as connecting point)

        density (float): Density of the system

        builder (Builder): Builder (One of `EMC` or `PSP`)

        natoms_total (int): Total number of atoms in the system, one of
                            this attribute and `nchains_total` has to be
                            provided but not both (providing both will result
                            in an error); default: `None`

        nchains_total (int): Total number of polymer chains in the system, one
                             of this attribute and `natoms_total` has to be
                             provided but not both (providing both will result
                             in an error); default: `None`

        natoms_per_chain (int): Number of atoms per polymer chain, one of this
                                attribute, `mw_per_chain`, and `ru_per_chain`
                                has to be provided but not more than 1
                                (providing more than 1 will result in an
                                error); default: `None`

        mw_per_chain (int): Molecular weight of the polymer, one of this
                            attribute, `natoms_per_chain`, and `ru_per_chain`
                            has to be provided but not more than 1
                            (providing more than 1 will result in an
                            error); default: `None`

        ru_per_chain (int): Number of repeating unit per polymer chain, one of
                            this attribute, `natoms_per_chain`, and
                            `mw_per_chain` has to be provided but not more
                            than 1 (providing more than 1 will result in an
                            error); default: `None`

        end_cap_smiles (str): SMILES string of the end-cap unit for polymers
                              ; default: `"*C"` (hint: put `"*[H]"` for end
                              capping with -H)

        data_fname (str): File name of the output data file, which will be
                          read in by LAMMPS
                          [read_data](https://docs.lammps.org/read_data.html)
                          command; default: `"data.lmps"`
    '''

    def __init__(self,
                 smiles: str,
                 density: float,
                 builder: Builder,
                 *,
                 natoms_total: Optional[int] = None,
                 nchains_total: Optional[int] = None,
                 natoms_per_chain: Optional[int] = None,
                 mw_per_chain: Optional[int] = None,
                 ru_per_chain: Optional[int] = None,
                 end_cap_smiles: str = '*C',
                 data_fname: str = 'data.lmps'):

        self._smiles = smiles
        self._density = density
        self._natoms_total = natoms_total
        self._nchains_total = nchains_total
        self._builder = builder
        self._mw_per_chain = mw_per_chain
        self._natoms_per_chain = natoms_per_chain
        self._ru_per_chain = ru_per_chain
        self._end_cap_smiles = end_cap_smiles
        self._data_fname = data_fname

        # Make sure only 1 system size option is given
        validate_options(self, SYSTEM_SIZE_OPTIONS)
        # Make sure only 1 chain length option is given
        validate_options(self, CHAIN_LENGTH_OPTIONS)
        # Make sure the provided builder is correct
        self._validate_builder()
        # Calculate system specs such as chain length, # of polymers
        self._calculate_system_spec()

    def __repr__(self) -> str:
        return type(self).__name__

    @property
    def smiles(self) -> str:
        return self._smiles

    @property
    def end_cap_smiles(self) -> str:
        return self._end_cap_smiles

    @property
    def data_fname(self) -> str:
        return self._data_fname

    @property
    def builder(self) -> Builder:
        return self._builder

    @smiles.setter
    def smiles(self, smiles: str):
        self._smiles = smiles
        self._calculate_system_spec()

    @end_cap_smiles.setter
    def end_cap_smiles(self, end_cap_smiles: str):
        self._end_cap_smiles = end_cap_smiles
        self._calculate_system_spec()

    @builder.setter
    def builder(self, builder: str):
        self._builder = builder
        self._validate_builder()

    def _validate_builder(self) -> None:
        if not isinstance(self._builder, Builder):
            raise ValueError('Invalid builder value, please provide '
                             'either a EMC or PSP object')

    def _calculate_polymer_spec(self) -> Tuple[int]:
        # Get the number of atoms of a repeating unit and determine the polymer
        # chain length
        mol = Chem.MolFromSmiles(self._smiles)
        natoms_per_RU = mol.GetNumAtoms(onlyExplicit=0) - 2
        if self._natoms_per_chain:
            self._final_ru_per_chain = round(self._natoms_per_chain /
                                             natoms_per_RU)
        elif self._mw_per_chain:
            mw_per_RU = Descriptors.ExactMolWt(mol)
            self._final_ru_per_chain = round(self._mw_per_chain / mw_per_RU)
        else:
            self._final_ru_per_chain = self._ru_per_chain

        # Get the number of atoms of a end-cap unit
        mol_end_cap = Chem.MolFromSmiles(self._end_cap_smiles)
        natoms_end_cap = mol_end_cap.GetNumAtoms(onlyExplicit=0) - 1
        return natoms_per_RU, natoms_end_cap

    def _calculate_system_spec(self) -> None:
        natoms_per_RU, natoms_end_cap = self._calculate_polymer_spec()

        # calculate final total number of chains
        if self._natoms_total:
            self._final_nchains_total = round(
                self._natoms_total /
                (natoms_per_RU * self._final_ru_per_chain +
                 natoms_end_cap * 2))
        else:
            self._final_nchains_total = self._nchains_total

        # calculate final total number of atoms
        self._final_natoms_total = (
            self._final_ru_per_chain * natoms_per_RU +
            natoms_end_cap * 2) * self._final_nchains_total

        Pmdlogging.info('System stats generated\n'
                        '--------Polymer Stats--------\n'
                        f'SMILES: {self._smiles}\n'
                        f'Natom_per_RU: {natoms_per_RU}\n'
                        f'length: {self._final_ru_per_chain}\n'
                        f'Nchains: {self._final_nchains_total}\n'
                        f'Total number of atoms: {self._final_natoms_total}\n'
                        '-----------------------------')

    def write_data(self, output_dir: str = '.', cleanup: bool = True) -> None:
        '''Method to make LAMMPS data file (which contains coordinates and force
        field parameters)

        Parameters:
        output_dir (str): Directory for the generated LAMMPS data file
                          ; default: `"."`

        cleanup (bool): Whether to clean up files other than the LAMMPS data
                        file PSP generated

        Returns:
            None
        '''

        self._builder.write_data(output_dir, self._smiles, self._density,
                                 self._final_natoms_total,
                                 self._final_ru_per_chain,
                                 self._final_nchains_total,
                                 self._end_cap_smiles, self._data_fname,
                                 cleanup)


class SolventSystem(System):
    '''Template object to contain System with solvents initialization settings

    Attributes:
        smiles (str): SMILES string of the polymer (use * as connecting point)

        solvent_smiles (str): SMILES string of the solvent

        ru_nsolvent_ratio (float): The ratio of total number of repeating units
                                   in the system and total number of solvent
                                   molecules

        density (float): Density of the system

        builder (Builder): Builder (One of `EMC` or `PSP`)

        natoms_total (int): Total number of atoms in the system, one of
                            this attribute and `nchains_total` has to be
                            provided but not both (providing both will result
                            in an error); default: `None`

        nchains_total (int): Total number of polymer chains in the system, one
                             of this attribute and `natoms_total` has to be
                             provided but not both (providing both will result
                             in an error); default: `None`

        natoms_per_chain (int): Number of atoms per polymer chain, one of this
                                attribute, `mw_per_chain`, and `ru_per_chain`
                                has to be provided but not more than 1
                                (providing more than 1 will result in an
                                error); default: `None`

        mw_per_chain (int): Molecular weight of the polymer, one of this
                            attribute, `natoms_per_chain`, and `ru_per_chain`
                            has to be provided but not more than 1
                            (providing more than 1 will result in an
                            error); default: `None`

        ru_per_chain (int): Number of repeating unit per polymer chain, one of
                            this attribute, `natoms_per_chain`, and
                            `mw_per_chain` has to be provided but not more
                            than 1 (providing more than 1 will result in an
                            error); default: `None`

        end_cap_smiles (str): SMILES string of the end-cap unit for polymers
                              ; default: `"*C"` (hint: put `"*[H]"` for end
                              capping with -H)

        data_fname (str): File name of the output data file, which will be
                          read in by LAMMPS
                          [read_data](https://docs.lammps.org/read_data.html)
                          command; default: `"data.lmps"`
    '''

    def __init__(self,
                 smiles: str,
                 solvent_smiles: str,
                 ru_nsolvent_ratio: float,
                 density: float,
                 builder: Builder,
                 *,
                 natoms_total: Optional[int] = None,
                 nchains_total: Optional[int] = None,
                 natoms_per_chain: Optional[int] = None,
                 mw_per_chain: Optional[int] = None,
                 ru_per_chain: Optional[int] = None,
                 end_cap_smiles: str = '*C',
                 data_fname: str = 'data.lmps'):

        self._solvent_smiles = solvent_smiles
        self._ru_nsolvent_ratio = ru_nsolvent_ratio

        super().__init__(smiles,
                         density,
                         builder,
                         natoms_total=natoms_total,
                         nchains_total=nchains_total,
                         natoms_per_chain=natoms_per_chain,
                         mw_per_chain=mw_per_chain,
                         ru_per_chain=ru_per_chain,
                         end_cap_smiles=end_cap_smiles,
                         data_fname=data_fname)

    @property
    def solvent_group(self):
        return f'molecule <= {self._nsolvents}'

    @property
    def polymer_group(self):
        return f'molecule > {self._nsolvents}'

    def _calculate_system_spec(self) -> None:
        natoms_per_RU, natoms_end_cap = self._calculate_polymer_spec()

        # Get the number of atoms of a solvent molecule
        mol_solvent = Chem.MolFromSmiles(self._solvent_smiles)
        natoms_solvent = mol_solvent.GetNumAtoms(onlyExplicit=0)

        # Calculate number of polymer chains and solvents based on target total
        # number of atoms
        natoms_total_onechain = (
            self._ru_nsolvent_ratio * self._final_ru_per_chain * natoms_solvent
        ) + (self._final_ru_per_chain * natoms_per_RU + natoms_end_cap * 2)

        # calculate final total number of chains
        if self._natoms_total:
            self._final_nchains_total = round(self._natoms_total /
                                              natoms_total_onechain)
        else:
            self._final_nchains_total = self._nchains_total

        # calculate final total number of solvents
        self._nsolvents = round(self._ru_nsolvent_ratio *
                                self._final_ru_per_chain *
                                self._final_nchains_total)

        # calculate final total number of atoms
        self._final_natoms_total = self._nsolvents * natoms_solvent + (
            self._final_ru_per_chain * natoms_per_RU +
            natoms_end_cap * 2) * self._final_nchains_total

        # Calculate extra stats for logging use
        final_nsol_nRU_ratio = self._nsolvents / (self._final_ru_per_chain *
                                                  self._final_nchains_total)

        Pmdlogging.info(
            'System stats generated\n'
            '--------Polymer Stats--------\n'
            f'Polymer SMILES: {self._smiles}\n'
            f'Polymer length: {self._final_ru_per_chain}\n'
            f'Polymer Nchains: {self._final_nchains_total}\n\n'
            '--------Solvent Stats--------\n'
            f'Solvent SMILES: {self._solvent_smiles}\n'
            f'Solvent number: {self._nsolvents}\n\n'
            '--------System Stats---------\n'
            f'Target Nsolvents/Nrepeatunits: {self._ru_nsolvent_ratio}\n'
            f'Final Nsolvents/Nrepeatunits: {final_nsol_nRU_ratio}\n'
            f'Total number of atoms: {self._natoms_total}\n'
            '-----------------------------')

    def write_data(self, output_dir: str = '.', cleanup: bool = True) -> None:
        '''Method to make LAMMPS data file (which contains coordinates and force
        field parameters)

        Parameters:
        output_dir (str): Directory for the generated LAMMPS data file
                          ; default: `"."`

        cleanup (bool): Whether to clean up files other than the LAMMPS data
                        file PSP generated

        Returns:
            None
        '''

        self._builder.write_solvent_data(
            output_dir, self._smiles, self._solvent_smiles, self._density,
            self._final_natoms_total, self._final_ru_per_chain,
            self._nsolvents, self._final_nchains_total, self._end_cap_smiles,
            self._data_fname, cleanup)
