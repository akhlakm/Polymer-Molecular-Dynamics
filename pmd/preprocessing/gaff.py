"""
Create a polymer system with GAFF2 as force field using the following steps:
- Build a small 1 chain system with OPLS-AA force field using EMC.
- Calculate AM1-BCC and GAFF2 atom types of the chain with AnteChamber.
- Create a mapping between the OPLS and GAFF2 atom types with OpenMol.
- Build the target larger polymer system using EMC.
- Using the mapping, convert the atom types and charges to GAFF2 types.
- Build the polymer system using tLeap for Amber.
- Convert the Amber system into LAMMPS data using OpenMol.

Requirements:
- emc-pypi (pip install emc-pypi)
- Ambertools (conda install -c conda-forge ambertools)
- OpenMol  (pip install openmol >= 1.2.1)
   
Author: Akhlak Mahmood
Version: 2023-07-30

"""

import os
import shutil

import openmol.tripos_mol2
import openmol.amber_parm7
import openmol.lammps_full

from pmd.util import Pmdlogging
from pmd.preprocessing.emc import EMCTool
from pmd.preprocessing.amber import AmberTool
from pmd.preprocessing.maps import EMC2GAFF



class GAFF2:
    """ Build a polymer system using GAFF2 as a force field.
    
    Attributes:
        working_dir (str)   Output directory of the files.
        data_fname  (str)   Lammps data file name.
    """
    def __init__(self, working_dir : str = ".", data_fname : str = "data.lmps"):
        self.emc_ff = "opls-aa"
        self.work_dir = working_dir
        self.temp_dir = os.path.join(self.work_dir, "_gaff")
        self.mapping = EMC2GAFF(self.emc_ff, self.temp_dir)
        self.map_file = os.path.join(self.work_dir, self.emc_ff+"2gaff.map")
        self.lammps_data = os.path.join(self.work_dir, data_fname)


    def load_or_create_mapping(self, smiles, end_cap_smiles,
                               cleanup : bool = False):
        """ Load or create a new force field mapping file for GAFF2. """
        if os.path.isfile(self.map_file):
            self.mapping.load_mapping(self.map_file)
        else:
            emc_prefix = "1chain"
            os.makedirs(self.temp_dir, exist_ok=True)

            # Create a small polymer chain with EMC.
            emc = EMCTool(self.emc_ff)
            emc.create_polymer_system(self.temp_dir, emc_prefix, smiles,
                                      1.0, end_cap_smiles, 100, 10, 1)
            emc.extract_files(self.temp_dir)

            # convert to mol2 to store atomic charges.
            self._emc_to_mol2(emc_prefix)
            self._calc_gaff2_charges(emc_prefix)

            # create a mapping with GAFF2
            emc_psf = os.path.join(self.temp_dir, emc_prefix+".psf")
            gaff_mol2 = os.path.join(self.temp_dir, emc_prefix+".gaff2.mol2")
            self.mapping.create(emc_psf, gaff_mol2, self.map_file)

            if cleanup:
                shutil.rmtree(self.temp_dir, ignore_errors=True)


    def create_system(self, smiles, density, natoms_total, length,
                      nchains, end_cap_smiles, cleanup=True):
        """ Create polymer system with GAFF2 using a FF mapping. """

        os.makedirs(self.temp_dir, exist_ok=True)
        emc_prefix = "system"

        emc = EMCTool(self.emc_ff)
        emc.create_polymer_system(self.temp_dir, emc_prefix, smiles, density,
                                  end_cap_smiles, natoms_total, length, nchains)
        emc.extract_files(self.temp_dir)

        self._emc_to_mol2(emc_prefix)
        self._remap_emc_to_gaff(emc_prefix)
        self._build_amber_system(emc_prefix)
        self._write_lammps_data(emc_prefix)
        self._copy_important_files(emc_prefix)

        if cleanup:
            shutil.rmtree(self.temp_dir, ignore_errors=True)


    def _emc_to_mol2(self, emc_prefix):
        """ Convert the EMC generated prefix.pdb to prefix.mol2. """
        input_pdb = emc_prefix+".pdb"
        input_psf = emc_prefix+".psf"
        out_mol2 =  emc_prefix+".mol2"

        amber = AmberTool(self.temp_dir)
        amber.convert_format(input_pdb, out_mol2, topology=input_psf)

        return out_mol2


    def _calc_gaff2_charges(self, emc_prefix):
        """ Calculate GAFF specific atom types and charges using AnteChamber."""
        amber = AmberTool(output_dir=self.temp_dir)
        amber.calc_atomtypes_charges(emc_prefix+".mol2",
                                     emc_prefix+".gaff2.mol2", resname="POL",
                                     quiet=True, cleanup=False)


    def _remap_emc_to_gaff(self, emc_prefix):
        """ Update the atom and charges of EMC built system mol2 file with
        GAFF2 parameters. """
        
        if self.mapping is None:
            raise ValueError("mapping not loaded/created")
        
        input_mol2 = os.path.join(self.temp_dir,
                                  emc_prefix+".mol2")
        out_mol2 = os.path.join(self.temp_dir, emc_prefix+".gaff2.mol2")
        system = openmol.tripos_mol2.read(input_mol2)

        # Set atom type and charge
        for i, orig_atom_name in enumerate(system.atom_name):
            system.atom_type[i] = self.mapping.ff2_name2type[orig_atom_name]
            system.atom_q[i] = self.mapping.ff2_name2charge[orig_atom_name]

        system = openmol.tripos_mol2.build(system)
        writer = openmol.tripos_mol2.Writer(system, out_mol2)
        writer.write()
        return out_mol2


    def _build_amber_system(self, emc_prefix, water="tip3p"):
        """ Build Amber simulation files from a MOL2 file. """
        input_leap = emc_prefix+'.gaff2.tleap'
        input_mol2 = emc_prefix+".gaff2.mol2"
        out_prmtop = emc_prefix+".gaff2.prmtop"
        out_restrt = emc_prefix+".gaff2.rst7"

        with open(os.path.join(self.temp_dir, input_leap), 'w+') as fp:
            fp.write(f"# tleap script to build polymer system.\n")
            fp.write(f"source leaprc.gaff2\n")
            fp.write(f"source leaprc.water.{water}\n")
            fp.write(f"logFile {emc_prefix}.leap.log\n")
            fp.write(f"sys = loadMol2 {input_mol2}\n")
            fp.write(f"addIons2 sys Na+ 0\n")
            fp.write(f"addIons2 sys Cl- 0\n")
            fp.write(f"setBox sys vdw 1\n")
            fp.write(f"\n")
            fp.write(f"saveAmberParm sys {out_prmtop} {out_restrt}\n")
            fp.write(f"charge sys\n")
            fp.write(f"quit")

        # Note!! EMC does not wrap the polymer chains.
        # This will create a huge box. We will manually set boxsize
        # during coversion to LAMMPS data.

        amber = AmberTool(self.temp_dir)
        amber.build_amber_system(input_leap, quiet=True)

        if os.path.isfile(out_prmtop):
            Pmdlogging.info(f"Preprocess - {out_prmtop} generated")

        if os.path.isfile(out_restrt):
            Pmdlogging.info(f"Preprocess - {out_restrt} generated")


    def _update_box_from_emc(self, emc_prefix : str, mol : dict):
        """ Update the box information of the OpenMol object from EMC PDB. """
        emc_pdb = os.path.join(self.temp_dir, emc_prefix+".pdb")
        # Read the emc pdb file to parse boxsize.
        # CRYST1  41.2179  41.2179  41.2179     90     90     90 P 1          1
        with open(emc_pdb) as fp:
            for line in fp:
                if "CRYST1" in line:
                    break

        if "CRYST1" not in line:
            raise ValueError("No box information found in EMC PDB.")
        
        words = line.split()
        mol['box_x'] = float(words[1])
        mol['box_y'] = float(words[2])
        mol['box_z'] = float(words[2])
        mol['box_alpha'] = float(words[4])
        mol['box_beta'] = float(words[5])
        mol['box_gamma'] = float(words[6])

        return mol


    def _write_lammps_data(self, emc_prefix):
        """ Convert amber simulation files into lammps data using OpenMol. """
        prmtop = os.path.join(self.temp_dir, emc_prefix+".gaff2.prmtop")
        restrt = os.path.join(self.temp_dir, emc_prefix+".gaff2.rst7")

        # read amber parm files
        p = openmol.amber_parm7.read(prmtop, restrt)

        # calculate necessary amber items
        p = openmol.amber_parm7.build(p)

        # calculate necessary lammps items
        p = openmol.lammps_full.build(p)

        # set title
        p.title = f"Polymer System with GAFF2"

        # use original emc box info
        p = self._update_box_from_emc(emc_prefix, p)

        # write lammps data file
        lmp = openmol.lammps_full.Writer(p, self.lammps_data)
        lmp.write()


    def _copy_important_files(self, emc_prefix):
        prmtop = os.path.join(self.temp_dir, emc_prefix+".gaff2.prmtop")
        prmtop_to = os.path.join(self.work_dir, emc_prefix+".prmtop")
        shutil.copy2(prmtop, prmtop_to)


if __name__ == "__main__":
    # Test environment, if executed directly.

    gaff = GAFF2(working_dir="tools_tests")
    gaff.load_or_create_mapping("[*]CC[*]", "*C", cleanup=True)
    gaff.create_system(smiles='*CC*', density=1.0, end_cap_smiles="*C",
                       natoms_total=1000, length=50, nchains=20)
