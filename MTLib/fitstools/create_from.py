
from astropy.io import fits

class create_from:
    '''
    Create a new fits file with path+name <to_file> from the <from_file> using a with operation.
    The hdu is returned from the with operation, so data can be modified before creating the new file.
    '''

    def __init__(self, from_file: str, to_file:str, index:int=0):
        self.from_file = from_file
        self.to_file = to_file
        self.index = index
        self.hdu=None

    def __enter__(self):
        '''Load the HDU and the WCS'''
        self.hdu = fits.open(self.from_file)[self.index]
        return self.hdu

    def __exit__(self, exc_type, exc_value, traceback):
        self.hdu.writeto(self.to_file, overwrite=True)