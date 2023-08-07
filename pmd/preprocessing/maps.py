import os
import OpenMol
import OpenMol.psf
import OpenMol.mapping
from pmd.util import Pmdlogging

class EMC2GAFF(OpenMol.mapping.FFMap):
    """ Mapping between EMC generated OPLSAA and GAFF2. """
    def __init__(self, emc_ff : str, output_dir : str = ".") -> None:
        super().__init__(emc_ff, "gaff2")
        self.work_dir = output_dir
        os.makedirs(self.work_dir, exist_ok=True)


    def _read_emc_psf(self, input_psf : str, output : str):
        """ Read the EMC generated system.psf.gz and save as json. """
        psf_file = OpenMol.psf.Reader()
        psf_file.read(input_psf)
        OpenMol.write_json(psf_file.Mol, output)
        Pmdlogging.info(f"Preprocess - {output} generated")


    def _create_ff_mapping(self, emc_psf_json : str, gaff_mol2_file : str,
                           map_file : str):
        """ Create mapping between amber and EMC/opls force field. """
        self.load_ff1_file(emc_psf_json)
        self.load_ff2_file(gaff_mol2_file)
        self.save_mapping(map_file)


    def create(self, emc_psf : str, gaff_mol2 : str, map_file : str):
        psf_json = emc_psf[:-4] + ".json"
        self._read_emc_psf(emc_psf, psf_json)
        self._create_ff_mapping(psf_json, gaff_mol2, map_file)
