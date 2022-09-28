import pmd

if __name__ == '__main__':
    # Build a Polyethlyene system
    system = pmd.System(smiles='*CC*',
                        density=0.8,
                        natoms_total=10000,
                        natoms_per_chain=600,
                        builder=pmd.EMC(force_field='pcff'))

    # Equilibration + HeatFluxMeasurement
    # (will print thermal conductivity result at the end)
    lmp = pmd.Lammps(read_data_from=system)
    lmp.add_procedure(pmd.Minimization())
    lmp.add_procedure(
        pmd.Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163))
    lmp.add_procedure(
        pmd.HeatFluxMeasurement(duration=10**7,
                                T=300,
                                dump_image=True,
                                reset_timestep_before_run=True))

    # Setup for the Torque scheduling system's job file
    job = pmd.Torque(run_lammps=lmp,
                     jobname='PE_heatflux_measurement',
                     project='GT-rramprasad3-CODA20',
                     nodes=2,
                     ppn=24,
                     walltime='24:00:00')

    # Create all the files in the PE_heatflux_measurement folder
    run = pmd.Pmd(system=system, lammps=lmp, job=job)
    run.create(output_dir='PE_heatflux_measurement', save_config=True)
