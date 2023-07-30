import pmd

if __name__ == '__main__':
    # Build a Polyethlyene system
    system = pmd.System(smiles='*CC*',
                        density=0.8,
                        natoms_total=10000,
                        natoms_per_chain=600,
                        builder=pmd.OpenMol(force_field='gaff2-am1bcc'))

    # Equilibration + Uniaxial tensile deformation
    # We specify the dcd output format which can be read by VMD
    # and CppTraj for further postprocessing.
    lmp = pmd.Lammps(read_data_from=system, use_gpu=True)
    lmp.add_procedure(pmd.Minimization(maxiter=5000, maxeval=10000,
                                       dump_every=100,
                                       dump_fname="min.dcd"))
    lmp.add_procedure(
        pmd.Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163,
                          dump_fname="equil.dcd"))
    lmp.add_procedure(
        pmd.TensileDeformation(
            duration=10**7,  # [fs]
            erate=10**-6,  # [fs]
            T=300,
            P=1,
            dump_fname="prod.dcd",
            reset_timestep_before_run=True))


    # Create all the files in the current folder
    run = pmd.Pmd(system=system, lammps=lmp)
    run.create(output_dir='.', save_config=True, cleanup=False)
