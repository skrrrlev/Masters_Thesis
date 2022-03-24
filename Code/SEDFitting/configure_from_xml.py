# imports
from astropy import units as U
import numpy as np
from os.path import isfile, isdir
import sys
from typing import NoReturn

# private imports
from C4S import Cataloguer
from MTLib.data.xml import Parser as xmlpar

def usage() -> NoReturn:
    docs = """
    python configure_from_xml.py <xml_file> [<output_path> [<output_name>]]
        <xml_file>    : path to xml file
            
        <output_path> : path to output
            If output_path is not given, the path will defeault to directory of the xml file.
        <output_name> : name of output
            If output_name is not given, the name will default to the name of the xml file.
    """
    print(docs)
    exit()

def setup() -> "tuple[str,str,str]":
    number_of_args = len(sys.argv) - 1
    if not number_of_args or number_of_args > 3:
        usage()
    print(number_of_args)

    xml_file: str = sys.argv[1].replace('\\','/')
    if not isfile(xml_file):
        usage()
    print(xml_file)

    if number_of_args > 1:
        output_path: str = sys.argv[2]
        if not isdir(output_path):
            print(f'{output_path} is not a directory')
            usage()
    else:
        output_path = '/'.join(xml_file.split('/')[:-1])
    print(output_path)

    if number_of_args == 3:
        output_name: str = sys.argv[3]
    else:
        output_name: str = xml_file.split('/')[-1].strip('.xml')
    print(output_name)
    
    return xml_file, output_path, output_name


def main():
    xml_file, output_path, output_name = setup()
    
    # Parse xml file
    targets = xmlpar.load_targets(xml_file)
    observations = xmlpar.load_observations(xml_file)

    # Create cataloguer instance
    cat = Cataloguer(output_name, output_path, U.mJy)

    # Create targets in catalogue
    for target in targets:
        cat.create_target(
            id = target['id'],
            ra = target['RA'],
            dec = target['DEC'],
            z = target['z']
        )

    # Add observations of targets
    for observation in observations:
        if observation['unit'] == 'ABmag':
            unit = 'Jy'
            mag = observation['flux']
            mag_err = observation['flux']

            flux = (mag * U.ABmag).to(U.Jy).value
            _df_dm = - ( 3631 * np.log(10) * np.exp(-mag * np.log(10) / 2.5) ) / 2.5
            error = np.sqrt((_df_dm * mag_err)**2)
        else:
            unit = observation['unit']
            flux = observation['flux']
            error = observation['error']
        
        if observation['type'] == 'F':
            cat.add_observation(
                id = int(observation['target']),
                name = observation['name'],
                f = flux,
                Δf = error,
                unit = unit,
                code = int(observation['typeval'])
            )
        elif observation['type'] == 'E':
            cat.add_observation(
                id = int(observation['target']),
                name = observation['name'],
                f = flux,
                Δf = error,
                unit = unit,
                λ = observation['typeval']
            )
        else:
            raise ValueError('Observation type was not F ("Filter") or E ("Extra band").')

    cat.save()



if __name__ == '__main__':
    main()