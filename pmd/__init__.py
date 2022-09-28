import sys

from .core import (EMC, NPT, NVT, PSP, Equilibration, HeatFluxMeasurement,
                   Lammps, Minimization, MSDMeasurement, Pmd, RgMeasurement,
                   ShearDeformation, Slurm, SolventSystem, System,
                   TensileDeformation, TgMeasurement, Torque)
from .postprocessing import (calculate_diffusivity, calculate_Tg,
                             read_lammpstrj, read_lammpstrj_by_type)

if sys.version_info[:2] >= (3, 8):
    from importlib.metadata import PackageNotFoundError, version
else:
    from importlib_metadata import PackageNotFoundError, version

try:
    dist_name = __name__
    __version__ = version(dist_name)
except PackageNotFoundError:
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError

__all__ = [
    Lammps, Pmd, Torque, Slurm, System, SolventSystem, EMC, PSP, Minimization,
    Equilibration, RgMeasurement, TgMeasurement, MSDMeasurement,
    TensileDeformation, ShearDeformation, HeatFluxMeasurement, NVT, NPT,
    calculate_diffusivity, calculate_Tg, read_lammpstrj, read_lammpstrj_by_type
]
