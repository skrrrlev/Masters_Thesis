
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
        self.hdul=None
        self.hdu=None

    def __enter__(self):
        '''Load the HDU and the WCS'''
        self.hdul: fits.HDUList = fits.open(self.from_file)
        self.hdu: fits.PrimaryHDU = self.hdul[self.index]
        return self.hdu

    def __exit__(self, exc_type, exc_value, traceback):
        self.hdu.writeto(self.to_file, overwrite=True)
        self.hdul.close()

if __name__ == '__main__':
    with create_from(from_file='Output/Maps/s14_F140w/s14_F140w.fits', to_file='Output/Maps/s14_F140w/galaxy/s14_F140w.fits', index=0):
        pass