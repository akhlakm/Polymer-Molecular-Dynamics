import argparse
from typing import Optional, Sequence

from pmd.postprocessing import calculate_diffusivity, calculate_Tg
from pmd.util import Pmdlogging


def main(args: Optional[Sequence[str]] = None):
    args = parse_command_line(args)

    # TODO: setup a nice dict with all property options and
    # corresponding functions

    # For Tg calculation
    if args.property == 'Tg':
        if not args.result:
            Pmdlogging.error(
                'Please provide the Tg result file from the TgMeasurement '
                'procedure as the input argument, for example:'
                '\n$ pmd-analyze -p Tg temp_vs_density.txt')
            raise ValueError('Missing result file input argument')
        property_result = calculate_Tg(args.result)

    # For diffusivity calculation
    elif args.property == 'D':
        if args.result:
            property_result = calculate_diffusivity(args.result)
        else:
            property_result = calculate_diffusivity()

    return property_result


def parse_command_line(args: Optional[Sequence[str]] = None):
    parser = argparse.ArgumentParser(description=(
        'Analyze simulation result files using the Pmd Analysis module.'))
    parser.add_argument('result', help='Result file or folder', nargs='?')
    parser.add_argument('-p', '--property', type=str, help='Property name')
    args = parser.parse_args(args=args)

    return args


if __name__ == '__main__':
    main()
