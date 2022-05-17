#imports
from astropy import units as U
from os import listdir
from re import findall
import sys

from C4S import Cataloguer
from C4S.dataclasses import ColumnType

# private imports
from MTLib.data.stardust import CustomStardustFilters as csf
from MTLib.data.cosmos import SpecFile
from MTLib.data.cosmos import super_deblended as sd

cosmos_dir = 'Data/COSMOS/'
spec_input_dir = cosmos_dir+'COSMOS_uv_nir/'


def setup():
    name='opt_nir'
    path='Data/stardust/COSMOS/'
    
    number_of_input_args = len(sys.argv)-1
    if number_of_input_args == 1:
        name = sys.argv[1]
    elif number_of_input_args == 2:
        name = sys.argv[1]
        path = sys.argv[2]
    
    return name,path


def main():
    """ Main script """
    
    """Add custom filters to the filter set of Stardust. 
    They have to be updated if the stardust package has been reinstalled."""
    csf.add_file('Data/stardust/extra_filters.txt')
    
    name, path = setup()
    cat = Cataloguer(name,path)

    """Load .spec files"""
    # Load .spec files and add them to fits file structure
    specfiles = [SpecFile(spec_input_dir + file, unit=U.mJy) for file in listdir(spec_input_dir) if file.endswith('.spec')]
        

    """Load RA and DEC and create target dictionary"""
     # open .txt file containing redshifts and RADEC, and read the contents into a string
    with open(cosmos_dir+'coord_cosmos_id2020.txt','r') as f:
        content = f.read()

    # Add them to lists.
    targets = {}
    for specfile in specfiles:
        targets[specfile.id] = {}

        ra_pattern = r'(\d+.\d+)\s+\d+.\d+\s+'+f'{specfile.id}'
        targets[specfile.id]['ra'] = float(findall(ra_pattern,content)[0])

        dec_pattern = r'\d+.\d+\s+(\d+.\d+)\s+'+f'{specfile.id}'
        targets[specfile.id]['dec'] = float(findall(dec_pattern,content)[0])

        targets[specfile.id]['z'] = specfile.z
        
    """Define targets in catalogue"""
    for target in targets:
        cat.create_target(id=target, ra=targets[target]['ra'], dec=targets[target]['dec'], z=targets[target]['z'])

    """Add observations from .spec files to catalogue"""
    for specfile in specfiles:
        for observation in specfile.observations:
            cat.add_observation(
                id=specfile.id,
                name=observation.column,
                f=observation.flux,
                unit='mJy',
                Δf=observation.flux_error,
                code=observation.filter.stardust_code
            )
    '''
    """ Load deblended data from fits file """
    ir_data = sd.get_data(unit=U.mJy)
    for observation in ir_data:
        if observation.type == ColumnType.FILTER:
            cat.add_observation(
                id=observation.target_id,
                name=observation.column,
                f=observation.flux,
                Δf=observation.flux_error,
                unit='mJy',
                code=observation.filter.stardust_code
            )
        elif observation.type == ColumnType.EXTRA:
            cat.add_observation(
                id=observation.target_id,
                name=observation.column,
                f=observation.flux,
                Δf=observation.flux_error,
                unit='mJy',
                λ=observation.wavelength
            )
        else:
            raise ValueError('The type of the observation was not well defined')

    '''

    """ Save stardust components """
    cat.save()


if __name__ == '__main__':
    main()

    
    


