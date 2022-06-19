import inquirer


def main():
    questions = [
        inquirer.List(
            'system',
            message="What system do you need?",
            choices=[
                '1. System (amorphous homopolymer)',
                '2. SolventSystem (homopolymer + solvent)',
                '3. GasSystem (homopolymer + gas)'
            ],
        ),
        inquirer.List(
            'system_size',
            message="How do you want to determine system size?",
            choices=[
                '1. By total number of atoms',
                '2. By total number of polymer chains'
            ],
        ),
        inquirer.List(
            'chain_length',
            message="How do you want to determine polymer chain length?",
            choices=[
                '1. By number of atoms per chain',
                '2. By number of repeating units per chain',
                '3. By polymer molecular weight'
            ],
        ),
        inquirer.List(
            'builder',
            message=
            "What force field (Builder) do you want to use for this system?",
            choices=[
                '1. opls-lbcc (PSP)', '2. opls-cm1a (PSP)',
                '3. gaff2-gasteiger (PSP)', '4. gaff2-am1bcc (PSP)',
                '5. pcff (EMC)', '6. opls-aa (EMC)', '7. opls-ua (EMC)',
                '8. trappe (EMC)'
            ],
        ),
        inquirer.List(
            'lammps',
            message="What property do you want to compute?",
            choices=[
                '1. Glass transition temperature',
                '2. Gas/solvent diffusivity', '3. Viscosity',
                '4. Mechanical properties', '5. Thermal conductivity'
            ],
        ),
        inquirer.List(
            'job',
            message="What job scheduling system do you use?",
            choices=['1. Torque', '2. Slurm', '3. N/A (run locally)'],
        ),
    ]
    answers = inquirer.prompt(questions)
    print(answers)


if __name__ == '__main__':
    main()
