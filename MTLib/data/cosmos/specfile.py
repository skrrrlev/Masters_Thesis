# imports
import numpy as np
from astropy import units as U
from astropy.units import Quantity

# relative imports
from ..stardust.map_stardust_dotspec import get_cosmos_dotspec_map
from ..dataclasses import Observation

class SpecFile:
    '''
    Contains class methods for reading values from the .spec files,
    but can also be initialised with a .spec file, to automatically fill in all values.
    '''
    cs_map = get_cosmos_dotspec_map()
    _INVALID_ = -999

    def __init__(self, file: str, unit: Quantity =U.mJy) -> None:
        self.file = file
        self.id, self.z = self.get_header(file)
        self.observations = self.get_observations(self.file, unit)

    @classmethod
    def get_AB_mag(cls, file: str, start: int=13, end: int=49) -> "tuple[np.ndarray,np.ndarray]":
        '''Input a .spec file and get the magnitude and magnitude error columns.'''

        # open file and read it
        with open(file, 'r') as f:
            content = f.read()
        # generate list of lines in the file
        line_list = content.split('\n')
        # Get the values in each line of the lines list.
        value_array = np.array([[float(value) for value in line.split(' ') if value] for i,line in enumerate(line_list) if i>=start and i<=end])

        #extract magnitude and uncertainty from array
        mag = value_array[:,0] * U.ABmag
        mag_err = value_array[:,1] * U.ABmag
        return mag, mag_err

    @classmethod
    def __get_flux_uncertainty(cls, mag, mag_err):
        '''Get uncertainty in Jy from magnitude and magnitude uncertainty.'''
        df_dm = - ( 3631 * np.log(10) * np.exp(-mag * np.log(10) / 2.5) ) / 2.5
        return np.sqrt((df_dm * mag_err)**2) * U.Jy


    @classmethod
    def get_flux(cls, file:str, start: int=13, end: int=49, unit=U.mJy) -> "tuple[np.ndarray,np.ndarray]":
        '''Input a .spec file and get the flux and flux-error columns'''
        # get magnitude and uncertainty from specfile.
        # When data is missingm the magnitude field will be -999.
        mag,mag_err= SpecFile.get_AB_mag(file, start, end)
        # convert to specified unit
        # the -999 fields will become infinity (np.inf)
        flux = mag.to(unit)
        # get the flux uncertainty
        flux_err = SpecFile.__get_flux_uncertainty(mag.value,mag_err.value).to(unit)

        # treat the invalid values
        mask = (flux == np.inf)
        flux_err[mask] = SpecFile._INVALID_ * unit
        flux[mask] = SpecFile._INVALID_ * unit
        
        for i in range(len(flux)):
            if flux[i].value != SpecFile._INVALID_:
                print(f'{mag[i]:.2e} ± {mag_err[i]:.2e} --> {flux[i]:.2e} ± {flux_err[i]:.2e}\t SNR: {flux[i]/flux_err[i]:.3f}')
        
        return flux,flux_err

    @classmethod
    def get_observations(cls, file: str, unit: Quantity =U.mJy) -> "list[Observation]":
        '''Input a .spec file and get a list of all the observations within it.'''

        index = SpecFile.get_id(file)

        # extract the flux and the uncertainty on the flux
        flux,flux_err = SpecFile.get_flux(file, unit=unit)

        # Number of observations
        N = len(flux)
        
        # Add all the observations.
        obs: "list[Observation]" = []
        for i in range(N):
            if flux[i].value == SpecFile._INVALID_:
                continue
            
            new_obs = Observation(
                target_id = index,
                flux = flux[i].value,
                flux_error = flux_err[i].value,
                unit = unit,
                filter = SpecFile.cs_map[i]
            )
            obs.append(new_obs)
        
        return obs

    @classmethod
    def get_header(cls, file):
        '''Load the target id, spectroscopic redshift and photometric redshift from the header of the .spec file.'''
        with open(file, 'r') as f:
            line = f.read().split('\n')[1].split(' ')
        index = int(line[0])
        zspec = float(line[1])
        zphot = float(line[2])

        if zspec != -999:
            return index, zspec
        elif zphot != -999:
            return index,zphot
        else:
            return index,-99

    def get_zphot(self):
        '''Load the target id, spectroscopic redshift and photometric redshift from the header of the .spec file.'''
        with open(self.file, 'r') as f:
            line = f.read().split('\n')[1].split(' ')
        zphot = float(line[2])
        return zphot
    

    @classmethod
    def get_id(cls, file):
        index,_ = SpecFile.get_header(file)
        return index

    def get_target_id(self):
        return self.id