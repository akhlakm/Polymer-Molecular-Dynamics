import pmd

if __name__ == '__main__':
    # Build a system with a single 50-mer PE chain with 500 CO molecules
    system = pmd.SolventSystem(smiles='*CC*',
                               solvent_smiles='CO',
                               ru_nsolvent_ratio=float(500 / 50),
                               density=0.85,
                               nchains_total=1,
                               ru_per_chain=50,
                               builder=pmd.EMC(force_field='pcff'))

    # Equilibration + RgMeasurement
    # (will output Rg result at the end)
    lmp = pmd.Lammps(read_data_from=system)
    lmp.add_procedure(pmd.Minimization())
    lmp.add_procedure(
        pmd.NVT(duration=10**6,
                Tinit=300,
                Tfinal=300,
                reset_timestep_before_run=True))
    lmp.add_procedure(
        pmd.RgMeasurement(duration=2 * 10**7,
                          T=300,
                          P=1,
                          group=system.polymer_group,
                          dump_image=True,
                          reset_timestep_before_run=True))

    # Setup for the Torque scheduling system's job file
    job = pmd.Torque(run_lammps=lmp,
                     jobname='PE_in_CO_Rg_measurement',
                     project='GT-rramprasad3-CODA20',
                     nodes=2,
                     ppn=24,
                     walltime='24:00:00')

    # Create all the files in the PE_in_CO_Rg_measurement folder
    run = pmd.Pmd(system=system, lammps=lmp, job=job)
    run.create(output_dir='PE_in_CO_Rg_measurement', save_config=True)
