
import numpy as np
from astropy.io import fits
from ..files import extract_path, extract_filename
import astropy.units as u


def estimate_brightness(fits_file:str, sigma_file:str, mask_file:str) -> float:
    '''By specifying the fits file, the standard deviation file and the mask file, estimate the brightness of the source in ABmag.'''
    
    '''Load in the data from the files'''
    with fits.open(fits_file) as hdul:
        sci_data: np.ndarray = hdul[0].data
        photfnu = float(hdul[0].header['PHOTFNU'])
    with fits.open(sigma_file) as hdul:
        sigma_data: np.ndarray = hdul[1].data
    with fits.open(mask_file) as hdul:
        mask_data: np.ndarray = hdul[0].data
        source_mask = np.array(mask_data!=0,dtype=bool)
    
    '''Create a rough mask of the sources by taking all pixels with a std>=3'''
    rough_sources_pixel_mask = sci_data/sigma_data >= 3

    '''Extract all pixels that are in the rough sources pixel mask and is not in the mask of the other sources of the image'''
    rough_source_pixel_mask = np.logical_and(rough_sources_pixel_mask, np.logical_not(source_mask))
    pixels_to_sum = sci_data[rough_source_pixel_mask].flatten()

    '''Save the brightness mask so we can inspect it later.'''
    # Get the weight image in the fits file
    hdu = fits.open(fits_file)[0]
    hdu.header['EXTVER'] = 'BRIGHTNESS'

    # Set all pixels that are included to 1, and all that are not, to 0.
    hdu.data[rough_source_pixel_mask] = 1
    hdu.data[np.logical_not(rough_source_pixel_mask)] = 0

    # Now we save the brightness mask.
    output_folder = extract_path(fits_file).replace('Data','Output')
    output_name = extract_filename(fits_file) + '_brightness'
    output_file = output_folder + output_name + '.fits'
    print(f'Saving brightness_mask: {output_file}')
    hdu.writeto(output_file, overwrite=True)

    '''Calculate the estimate of the flux density and magnitude of the source'''
    flux_density = np.sum(pixels_to_sum)*photfnu * u.Jy
    magnitude = flux_density.to(u.ABmag)

    return magnitude.value


if __name__ == '__main__':
    print(estimate_brightness(
        fits_file='Output/Maps/HST/source14/galaxy/source14_F140w.fits',
        sigma_file='Output/Maps/HST/source14/galaxy/source14_F140w_sigma.fits',
        mask_file='Output/Maps/HST/source14/galaxy/source14_F140w_mask.fits'
    ))