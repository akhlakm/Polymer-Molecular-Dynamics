import os
import shutil

import openmol
import openmol.psf
import openmol.mapping
import openmol.tripos_mol2
import openmol.amber_parm7
import openmol.lammps_full


from pmd.core.Builder import EMC
from pmd.util import Pmdlogging
from pmd.preprocessing.utils import execute_command, extract_gzip


class GAFF2:
    def __init__(self, working_dir : str = ".", data_fname : str = "data.lmps"):
        self.mapping = None
        self.emc_ff = "opls-aa"
        self.work_dir = working_dir
        self.temp_dir = os.path.join(self.work_dir, "_gaff")
        self.map_file = os.path.join(self.work_dir, self.emc_ff+"2gaff.map")
        self.lammps_data = os.path.join(self.work_dir, data_fname)

    def load_or_create_mapping(self, smiles, end_cap_smiles,
                               cleanup : bool = False):
        """ Load or create a new force field mapping file for GAFF2. """
        if os.path.isfile(self.map_file):
            self.mapping = openmol.mapping.FFMap(self.emc_ff, "gaff2")
            self.mapping.load_mapping(self.map_file)

        else:
            emc_prefix = "1chain"
            self._check_antechamber()
            os.makedirs(self.temp_dir, exist_ok=True)

            self._create_chain(smiles, end_cap_smiles, emc_prefix)
            self._extract_emc()
            self._emc_to_mol2(emc_prefix)
            self._read_emc_psf(emc_prefix)
            self._calc_atomtypes_charges(emc_prefix)
            # self._generate_ff(emc_prefix)
            self._create_ff_mapping(emc_prefix)

            if cleanup:
                shutil.rmtree(self.temp_dir, ignore_errors=True)


    def create_system(self, smiles, density, natoms_total, length,
                      nchains, end_cap_smiles, cleanup=True):

        os.makedirs(self.temp_dir, exist_ok=True)
        emc_prefix = "system."

        # Generate the initial configuration of the system using EMC.
        emc = EMC(self.emc_ff)
        emc.write_data(self.temp_dir, smiles, density, natoms_total, length,
                       nchains, end_cap_smiles, f"{emc_prefix+self.emc_ff}.data",
                       cleanup=False)
        
        self._extract_emc()
        self._emc_to_mol2(emc_prefix+self.emc_ff)
        self._remap_emc_to_gaff(emc_prefix)
        self._build_amber_system(emc_prefix)
        self._write_lammps_data(emc_prefix)


    def _create_chain(self, smiles, end_cap_smiles, emc_prefix):
        """ Create a small polymer chain with EMC. """
        emc_data_file = emc_prefix + ".data"
        emc = EMC(self.emc_ff)
        emc.write_data(self.temp_dir, smiles, 1.0, 100, 10,
                       1, end_cap_smiles, emc_data_file, cleanup=False)

        # Make sure EMC ran successfully.
        if not os.path.isfile(os.path.join(self.temp_dir, emc_data_file)):
            raise RuntimeError("Failed to generate single chain.")

    def _extract_emc(self):
        input_pdb = os.path.join(self.temp_dir, "system.pdb.gz")
        input_psf = os.path.join(self.temp_dir, "system.psf.gz")
        if os.path.isfile(input_pdb):
            extract_gzip(input_pdb)
        if os.path.isfile(input_psf):
            extract_gzip(input_psf)


    def _emc_to_mol2(self, emc_prefix):
        """ Convert the EMC generated system.pdb.gz to mol2. """
        input_pdb = os.path.join(self.temp_dir, "system.pdb")
        input_psf = os.path.join(self.temp_dir, "system.psf")
        out_mol2 = os.path.join(self.temp_dir, emc_prefix+".mol2")

        cmd = f"cpptraj -p {input_psf} -y {input_pdb} -x {out_mol2}"

        Pmdlogging.info('Running cpptraj ...')
        res = execute_command(cmd, capture=False)

        if res.returncode != 0:
            raise RuntimeError("CppTraj failed.")

        if os.path.isfile(out_mol2):
            Pmdlogging.info(f"Preprocess - {out_mol2} generated")
        return out_mol2


    def _read_emc_psf(self, emc_prefix):
        """ Read the EMC generated system.psf.gz and save as json. """
        input_psf = os.path.join(self.temp_dir, "system.psf")
        out_json = os.path.join(self.temp_dir, emc_prefix+".json")

        # Extract
        if not os.path.isfile(input_psf):
            input_psf = extract_gzip(input_psf+".gz")

        psf_file = openmol.psf.Reader()
        psf_file.read(input_psf)
        openmol.write_json(psf_file.Mol, out_json)
        Pmdlogging.info(f"Preprocess - {out_json} generated")
        return out_json


    def _check_emc(self):
        #@todo: check if EMC is correctly installed.
        pass


    def _check_antechamber(self):
        antechamber = execute_command("antechamber -h").returncode == 0
        if not antechamber:
            raise RuntimeError("AnteChamber not installed. "\
                    "Please install ambertools using conda. "\
                    "conda install -c conda-forge ambertools ")


    def _check_tleap(self):
        tleap = execute_command("tleap -h").returncode == 0
        if not tleap:
            raise RuntimeError("tleap not installed. "\
                    "Please install ambertools using conda. "\
                    "conda install -c conda-forge ambertools ")


    def _calc_atomtypes_charges(self, emc_prefix, cleanup=False):
        """ Calculate GAFF specific atom types and charges using AnteChamber."""
        input_mol2 = emc_prefix+".mol2"
        out_mol2 = emc_prefix+".gaff2.mol2"

        atom_type = "gaff2"
        charge_type = "bcc"

        prev_wd = os.getcwd()
        os.chdir(self.temp_dir)

        cmd = f"antechamber -i {input_mol2} -fi mol2 " \
              f"-o {out_mol2} -fo mol2 " \
              f"-c {charge_type} -rn POL -at {atom_type} -s 0 -pl -1 " \
              f"-pf {'y' if cleanup else 'n'}"

        Pmdlogging.info('Running antechamber. This may take a while ...')
        res = execute_command(cmd, capture=False)

        if res.returncode != 0:
            raise RuntimeError("Antechamber failed.")

        if os.path.isfile(out_mol2):
            Pmdlogging.info(f"Preprocess - {out_mol2} generated")

        os.chdir(prev_wd)


    def _generate_ff(self, emc_prefix):
        """ Generate Amber force field file from a Mol2 file using parmchk2. """
        input_mol2 = emc_prefix+".gaff2.mol2"
        out_frcmod = emc_prefix+".gaff2.frcmod"

        atom_type = "gaff2"
        charge_type = "bcc"

        prev_wd = os.getcwd()
        os.chdir(self.temp_dir)

        cmd = f"parmchk2 -i {input_mol2} -o {out_frcmod} -f mol2 " \
              f"-s {2 if atom_type == 'gaff2' else 1 } -a Y -w Y "

        Pmdlogging.info('Running parmchk2 ...')
        res = execute_command(cmd, capture=False)

        if res.returncode != 0:
            raise RuntimeError("Parmchk2 failed.")

        if os.path.isfile(out_frcmod):
            Pmdlogging.info(f"Preprocess - {out_frcmod} generated")

        os.chdir(prev_wd)


    def _create_ff_mapping(self, emc_prefix):
        """ Create mapping between amber and EMC/opls force field. """

        psf_json = os.path.join(self.temp_dir, emc_prefix+".json")
        gaff_mol2 = os.path.join(self.temp_dir, emc_prefix+".gaff2.mol2")

        self.mapping = openmol.mapping.FFMap(self.emc_ff, "gaff2")
        self.mapping.load_ff1_file(psf_json)
        self.mapping.load_ff2_file(gaff_mol2)
        self.mapping.save_mapping(self.map_file)


    def _remap_emc_to_gaff(self, emc_prefix):
        """ Update the atom and charges of EMC built system mol2 file with
        GAFF2 parameters. """
        
        if self.mapping is None:
            raise ValueError("mapping not loaded/created")
        
        input_mol2 = os.path.join(self.temp_dir, emc_prefix+self.emc_ff+".mol2")
        out_mol2 = os.path.join(self.temp_dir, emc_prefix+"gaff2.mol2")
        system = openmol.tripos_mol2.read(input_mol2)

        # Set atom type and charge
        for i, orig_atom_name in enumerate(system.atom_name):
            system.atom_type[i] = self.mapping.ff2_name2type[orig_atom_name]
            system.atom_q[i] = self.mapping.ff2_name2charge[orig_atom_name]

        system = openmol.tripos_mol2.build(system)
        writer = openmol.tripos_mol2.Writer(system, out_mol2)
        writer.write()
        return out_mol2


    def _build_amber_system(self, emc_prefix):
        input_mol2 = emc_prefix+"gaff2.mol2"
        input_leap = emc_prefix+'gaff2.tleap'
        out_prmtop = emc_prefix+"gaff2.prmtop"
        out_restrt = emc_prefix+"gaff2.rst7"

        prev_wd = os.getcwd()
        os.chdir(self.temp_dir)

        with open(input_leap, 'w+') as fp:
            fp.write(f"# tleap script to build polymer system.\n")
            fp.write(f"logFile leap.log\n")
            fp.write(f"source leaprc.gaff2\n")
            fp.write(f"sys = loadMol2 {input_mol2}\n")
            fp.write(f"setBox sys centers 1\n")
            fp.write(f"saveAmberParm sys {out_prmtop} {out_restrt}\n")
            fp.write(f"quit")

        cmd = f"tleap -f {input_leap}"

        Pmdlogging.info('Running tleap. This may take a while ...')
        res = execute_command(cmd, capture=False)

        if res.returncode != 0:
            raise RuntimeError("tleap failed.")

        if os.path.isfile(out_prmtop):
            Pmdlogging.info(f"Preprocess - {out_prmtop} generated")

        if os.path.isfile(out_restrt):
            Pmdlogging.info(f"Preprocess - {out_restrt} generated")

        os.chdir(prev_wd)


    def _write_lammps_data(self, emc_prefix):
        prmtop = os.path.join(self.temp_dir, emc_prefix+"gaff2.prmtop")
        restrt = os.path.join(self.temp_dir, emc_prefix+"gaff2.rst7")

        # read amber parm files
        p = openmol.amber_parm7.read(prmtop, restrt)

        # calculate necessary amber items
        p = openmol.amber_parm7.build(p)

        # calculate necessary lammps items
        p = openmol.lammps_full.build(p)

        # set title
        p.title = f"Polymer System with GAFF2"

        # write lammps data file
        lmp = openmol.lammps_full.Writer(p, self.lammps_data)
        lmp.write()


if __name__ == "__main__":
    gaff = GAFF2(working_dir="gaff_tests")
    gaff.load_or_create_mapping("[*]CC[*]", "*C", cleanup=True)
    gaff.create_system(smiles='*CC*', density=1.0, end_cap_smiles="*C",
                       natoms_total=10000, length=500, nchains=20)
