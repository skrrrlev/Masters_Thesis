# imports
from astropy.units import Quantity
from astropy import units as U
from typing import Union

# relative imports
from .filter import Filter
from C4S.dataclasses import ColumnType

class Observation:
    '''Class containing a single observation.'''

    def __init__(self,target_id: int,flux: float,flux_error: float,unit: Quantity, filter: Union[Filter,None]=None,wavelength: Union[Quantity,None]=None):
        self.target_id = target_id
        self.flux = flux
        self.flux_error = flux_error
        self.unit = unit
        
        if isinstance(filter, Filter) and isinstance(wavelength, Quantity):
            raise ValueError('Cannot discern type of observation, since both filter and wavelength were specified.')
        elif isinstance(filter, Filter):
            self.type = ColumnType.FILTER
            self.filter = filter
            self.column = filter.get_generic_name()
        elif isinstance(wavelength, Quantity):
            self.type = ColumnType.EXTRA
            self.wavelength: float = wavelength.to(U.um).value
            self.column = f'{wavelength.value}{str(wavelength.unit)}'
        else:
            raise ValueError('Either specify filter or central wavelength, to create observation object')

    def set_name(self, name:str) -> None:
        self.column = name
    
    def get_name(self):
        return self.column