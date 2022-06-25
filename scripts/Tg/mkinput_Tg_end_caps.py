import pmd

if __name__ == '__main__':
    system = pmd.System(smiles='*COC*',
                        density=0.8,
                        natoms_total=5000,
                        natoms_per_chain=500,
                        builder=pmd.EMC(force_field='opls-aa'))

    lmp = pmd.Lammps(read_data_from=system)
    lmp.add_procedure(pmd.Minimization())
    lmp.add_procedure(
        pmd.Equilibration(Teq=600, Peq=1, Tmax=1000, Pmax=49346.163))
    lmp.add_procedure(
        pmd.TgMeasurement(Tinit=600,
                          Tfinal=100,
                          Tinterval=20,
                          step_duration=5000000))

    job = pmd.Torque(run_lammps=lmp,
                     jobname='test',
                     project='GT-rramprasad3-CODA20',
                     nodes=3,
                     ppn=24,
                     walltime='72:00:00')

    run = pmd.Pmd(system, lmp, job)

    try:
        # default end-cap SMILES is '*C'
        run.create('test', save_config=True)
    except Exception:
        system.end_cap_smiles = '*[H]'
        run.create('test', save_config=True)
