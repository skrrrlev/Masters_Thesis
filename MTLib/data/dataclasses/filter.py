
from dataclasses import dataclass, field
from astropy.units.quantity import Quantity

@dataclass(order=True,frozen=True)
class Filter:
    '''Dataclass describing a filter used for fitting purposes.'''
    instrument: str = field(compare=False)
    telescope: str = field(compare=False)
    band: str = field(compare=False)
    spec_wavelength: Quantity = field(compare=False, repr=False)
    spec_width: Quantity = field(compare=False, repr=False)
    stardust_code: int 

    def get_generic_name(self):
        return f'{self.telescope}_{self.instrument}_{self.band}'#/μ-{self.spec_wavelength.to(U.AA).value:.0f}Å/w-{self.spec_width.to(U.AA).value:.0f}Å'