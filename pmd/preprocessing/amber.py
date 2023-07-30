import os
from pmd.util import Pmdlogging
from pmd.preprocessing.utils import execute_command

class AmberTool:
    """ A low-level python interface to AmberTools. """
    def __init__(self, output_dir : str = ".") -> None:
        self._check_ambertools()
        os.makedirs(output_dir, exist_ok=True)
        self.work_dir = output_dir
        

    def _check_ambertools(self):
        antechamber = execute_command("antechamber -h").returncode == 0
        if not antechamber:
            raise RuntimeError("AnteChamber not found. "\
                    "Please install ambertools using conda. "\
                    "conda install -c conda-forge ambertools ")


    def _run_in_work_dir(self, command : str, save_prefix : str = None):
        capture = True if save_prefix else False        
        prev_wd = os.getcwd()

        try:
            os.chdir(self.work_dir)
            res = execute_command(command, capture=capture)
            if res.returncode != 0:
                raise RuntimeError(command)
            
            if save_prefix:
                with open(f"{save_prefix}.stdout.txt", "a+") as fp:
                    fp.write(res.stdout.decode("utf-8"))
                with open(f"{save_prefix}.stderr.txt", "a+") as fp:
                    fp.write(res.stderr.decode("utf-8"))

        finally:
            os.chdir(prev_wd)


    def calc_atomtypes_charges(self, input_mol2 : str, output_mol2 : str,
                               resname : str = "POL", quiet : bool = True,
                               netcharge : int = 0,
                               cleanup=False):
        """ Calculate GAFF specific atom types and charges using AnteChamber."""

        atom_type = "gaff2"
        charge_type = "bcc"

        cmd = f"antechamber -i {input_mol2} -fi mol2 " \
              f"-o {output_mol2} -fo mol2 " \
              f"-c {charge_type} -rn {resname} -at {atom_type} -s 0 -pl -1 " \
              f"-nc {netcharge} -pf {'y' if cleanup else 'n'}"

        Pmdlogging.info('Running antechamber. This may take a while ...')
        saveas = "antechamber" if quiet else None
        self._run_in_work_dir(cmd, saveas)


    def generate_ff_file(self, input_mol2 : str, output_ff_file : str,
                         quiet : bool = True):
        """ Generate Amber force field file from a Mol2 file using parmchk2. """
        atom_type = "gaff2"
        cmd = f"parmchk2 -i {input_mol2} -o {output_ff_file} -f mol2 " \
              f"-s {2 if atom_type == 'gaff2' else 1 } -a Y -w Y "

        saveas = "parmchk2" if quiet else None
        self._run_in_work_dir(cmd, saveas)


    def build_amber_system(self, input_script : str, quiet : bool = True):
        """ Build Amber simulation files using a TLeap script. """
        cmd = f"tleap -f {input_script}"
        Pmdlogging.info('Running tleap ...')
        saveas = "tleap" if quiet else None
        self._run_in_work_dir(cmd, saveas)


    def convert_format(self, trajectory : str, output_file : str,
                       topology : str = None, quiet : bool = True):
        """ Convert file formats using CppTraj. """
        if topology is None:
            cmd = f"cpptraj -y {trajectory} -x {output_file}"
        else:
            cmd = f"cpptraj -p {topology} -y {trajectory} -x {output_file}"

        Pmdlogging.info('Running cpptraj ...')
        saveas = "cpptraj.convert" if quiet else None
        self._run_in_work_dir(cmd, saveas)
