from typing import List, Dict, Any, Tuple, Union
from itertools import groupby
from csv import writer
from operator import itemgetter
from re import search, match
from os.path import exists

from atb_api import API

api = API(debug=True)

POSSIBLE_ANGLES = set([i * 10.0 for i in range(-36, 36 + 1)] + [None])

DATA_FILE_NAME = 'QM_MM_Gas_Phase_Torsion_Scan_Individual_Results_with_CCSD_T_CBS_baseline.sdf'

assert exists(DATA_FILE_NAME), '''ERROR:Can't find {0} in current directory'''.format(DATA_FILE_NAME)

def try_float(x: str) -> float:
    try:
        print(float(x))
        return float(x)
    except ValueError:
        return float('nan')

def molecule_name_and_dict_for(molecule_lines: List[str]) -> Dict[str, Any]:
    molecule_dict = {'atom_lines': [], 'bond_lines': []}

    for (i, line) in enumerate(molecule_lines):
        if i == 0:
            molecule_name = line
        elif i == 1:
            molecule_dict['id'] = line
        elif i == 2:
            molecule_dict['software'] = line
        else:
            if line.startswith('>'):
                molecule_dict[line.split()[-1].replace('<', '').replace('>', '')] = molecule_lines[i + 1].strip()
            elif len(line) in [70, 40]:
                molecule_dict['atom_lines'].append(line)
            elif len(line) == 22:
                molecule_dict['bond_lines'].append(line)
            else:
                if True:
                    pass
                else:
                    print(len(line), line)

    return molecule_dict

def on_angle(molecule: Dict[str, Any]) -> float:
    try:
        return float(molecule['ScanVar_1'])
    except TypeError:
        raise Exception(molecule)

def on_smiles(molecule: Dict[str, Any]) -> float:
    try:
        return molecule['SMILES']
    except TypeError:
        raise Exception(molecule)

def sdf_for_molecule(molecule: Dict[str, Any]) -> str:
    return ''.join(
        ['A\n' * 3]
        +
        molecule["atom_lines"]
        +
        molecule["bond_lines"]
        +
        ['M  END']
    )

if __name__ == '__main__':

    molecules = []

    molecule_lines = []
    next_line_is_molecule_name = True

    with open(DATA_FILE_NAME) as fh:
        for line in fh:
            if next_line_is_molecule_name:
                molecule_name = line.strip()
                molecule_lines.append(line.strip())
                next_line_is_molecule_name = False
            elif line.startswith('$$$$'):
                molecules.append(molecule_lines)
                molecule_lines = []
                next_line_is_molecule_name = True
                if False and len(molecules) == 2:
                    break
            else:
                molecule_lines.append(line)

    molecules = [
        molecule_name_and_dict_for(molecule_lines)
        for molecule_lines in molecules
    ]

    molecules_grouped_by_names = [
        (molecule_id, list(molecules))
        for (molecule_id, molecules) in groupby(
            sorted(molecules, key=on_smiles),
            key=on_smiles,
        )
    ]

    molecules_grouped_by_names_with_angles = [
        (
            molecule_id,
            len(molecules),
            len(
                {
                    angle: len(list(molecules_with_angle))
                    for (angle, molecules_with_angle) in groupby(
                        sorted(molecules, key=on_angle),
                        key=on_angle,
                    )
                },
            ),
        )
        for (molecule_id, molecules) in molecules_grouped_by_names

    ]

    with open('qm.csv', 'w') as fh:
        csv_writer = writer(fh)
        for row in sorted(molecules_grouped_by_names_with_angles, key=itemgetter(1), reverse=True):
            csv_writer.writerow(row)

    print({
        molecule_smiles: [
            match['molid']
            for match in
            api.Molecules.structure_search(structure_format='sdf', structure=sdf_for_molecule(molecules[0]), netcharge='*')['matches']
        ]
        for (molecule_smiles, molecules) in molecules_grouped_by_names
    })
