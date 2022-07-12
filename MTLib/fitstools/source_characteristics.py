
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, FK5, ICRS
from numpy import isnan
from typing import Union
from astropy import units as U

from matplotlib import pyplot as plt


def estimate_brightness(fits_file:str, sigma_file:str, mask_file:str, zero_point:float, brightness_mask_file:str = '') -> float:
    '''
    By specifying the fits file, the standard deviation file and the mask file, estimate the brightness of the source in mag.
    
    The <zero_point> parameter is the zero-point of the magnitude scale, as defined by the instrument and the scale.
    
    If the <brightness_mask_file> is not an empty string, this function will save the brightness mask so it can be inspected later on.

    Procedure:
    - Create a rough mask of the source pixels by selecting all pixels with a std>3
    - Create the mask by selecting all pixels in the rough mask, and all pixels not in the <mask_file>
    - Sum all pixels in the mask, and divide them by the 'EXPTIME' keyword from the map-header, to get the total photo-electrons/s.
    - Calculate the magnitude using m=-2.5*log10(flux-density) + <zero_point>
    '''
    
    '''Load in the data from the files'''
    with fits.open(fits_file) as hdul:
        sci_data: np.ndarray = hdul[0].data
        exptime = float(hdul[0].header['EXPTIME'])
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
    if brightness_mask_file:
        # Get the weight image in the fits file
        with fits.open(fits_file) as hdul:
            hdu = hdul[0]
            hdu.header['EXTVER'] = 'BRIGHTNESS'

            # Set all pixels that are included to 1, and all that are not, to 0.
            hdu.data[rough_source_pixel_mask] = 1
            hdu.data[np.logical_not(rough_source_pixel_mask)] = 0

            # Now we save the brightness mask.
            print(f'Saving brightness_mask: {brightness_mask_file}')
            hdu.writeto(brightness_mask_file, overwrite=True)

    '''Calculate the estimate of the flux density and magnitude of the source'''
    #flux_density = (np.sum(pixels_to_sum)*photfnu * U.Jy)#.to(U.Unit("erg cm^−2 s^−1 Hz^−1"))
    flux_density = np.sum(pixels_to_sum)/exptime
    magnitude = -2.5*np.log10(flux_density) + zero_point

    return magnitude

def get_pixel_in_map_from_coordinates(map:str, extension:int, ra:float, dec:float, frame=ICRS, unit='deg', anchor_quadrant:int=3) -> "tuple[float]":
    '''
    Get the zero-based pixel-index, from the edge of the <anchor_quadrant> of a map,
    that corresponds to the right ascension <ra> and the declination <dec>,
    using the <frame> in <unit>.
        - Default frame is astropy.coordinates.ICRS, and default unit is degrees ('deg')
    
    If coordinates are not within image bounds, a ValueError is raised.
    
    Returns (pixel_row, pixel_column).
    
    <anchor_quadrant> corresponds to pixel coordinates from the edge of the specified quadrant. 
    e.g. quadrant 1: top-right, quadrant 2: top-left, ...
    '''
    
    '''Load the HDU and the WCS'''
    hdul = fits.open(map)
    wcs = WCS(hdul[extension].header)
    rows,columns = hdul[extension].data.shape
    hdul.close()

    sky = SkyCoord(ra,dec,frame=frame, unit=unit)
    xp, yp = wcs.world_to_pixel(sky)
    if isnan(xp) or isnan(yp):
        raise ValueError('Input coordinates is not within image.')
    x = int(np.round(xp,0))
    y = int(np.round(yp,0))
    if anchor_quadrant == 1:
        return columns-x, rows-y
    elif anchor_quadrant == 2:
        return x,rows-y
    elif anchor_quadrant == 3:
        return x,y
    elif anchor_quadrant == 4:
        return columns-x, y
    else: 
        raise ValueError(f"Anchor quadrant = [1,2,3,4]. Recieved {anchor_quadrant}")


def get_pixel_in_map_from_brightest_pixel(map:str, extension:int, brightness_mask: str, anchor_quadrant:int=3) -> "tuple[float]":
    with fits.open(map) as hdul:
        data: np.ndarray = hdul[extension].data
        nrows, ncols = data.shape
        print(nrows, ncols)

    with fits.open(brightness_mask) as hdul:
        mask: np.ndarray = np.array(hdul[0].data,dtype='bool')

    masked_data = np.copy(data)
    masked_data[np.logical_not(mask)] = 0

    centre = np.unravel_index(np.argmax(masked_data),masked_data.shape) # row and column of brightest pixel
    bp_row: int = centre[0]+1 # traditionally the y-coordinate (not zero-based)
    bp_col: int = centre[1]+1 # traditionally the x-coordinate (not zero-based)
    print(f'brightest pixel from bottom-right: {bp_col}, {bp_row}')
    if anchor_quadrant == 1:
        return ncols+1-bp_col, nrows+1-bp_row
    elif anchor_quadrant == 2:
        return bp_col, nrows+1-bp_row
    elif anchor_quadrant == 3:
        return bp_col,bp_row
    elif anchor_quadrant == 4:
        return ncols+1-bp_col, bp_row 
    else: 
        raise ValueError(f"Anchor quadrant = [1,2,3,4]. Recieved {anchor_quadrant}")

def estimate_half_light_radius(map:str, extension:int, brightness_mask: str) -> float:
    '''
    Estimate the half-light radius of the source in pixels.
    Will only consider the pixels that are non-zero in the <brightness_mask>.
    
    Procedure:
    - Find the brightest pixel
    - Define an array of pixel distances: np.linspace(0,10,N) 
    - Define an array of angles: np.linspace(0,2π,M)
    - Loop through all of the angles
        - calculate the pixels corresponding to the different distances
        - get the value of the map at all of the distances
        - Add the values/M to an average-value-array
    - With the average-value-array, locate the index of the pixel that is closest the half-light radius
    - Extract the distance this corresponds to from the distance array.

    This method uses N=101, meaning that the minimum error that you should expect is ~0.1.

    This method uses M=25. A more precise approach is available if you know the P.A. of the source.
    '''

    with fits.open(map) as hdul:
        data: np.ndarray = hdul[extension].data

    with fits.open(brightness_mask) as hdul:
        mask: np.ndarray = np.array(hdul[0].data,dtype='bool')

    masked_data = np.copy(data)
    masked_data[np.logical_not(mask)] = 0

    centre = np.unravel_index(np.argmax(masked_data),masked_data.shape) # row and column of brightest pixel
    bp_row: int = centre[0]
    bp_col: int = centre[1]

    N = 101
    rs = np.linspace(0,10,N)
    ys = np.zeros(N)

    M = 25
    angles = np.linspace(0,2*np.pi,M)
    
    for angle in angles:
        
        def f(r:"Union[float,int,list,np.ndarray]") -> "Union[float,np.ndarray]":
            single_value_flag = False
            if isinstance(r,float) or isinstance(r,int):
                single_value_flag=True
                r = [r]
            r = np.array(r,dtype=float)

            col = np.cos(angle)*r+bp_col
            row = np.sin(angle)*r+bp_row

            col_indices = np.array(np.round(col,0),dtype=int)
            row_indices = np.array(np.round(row,0),dtype=int)

            values = masked_data[row_indices,col_indices]
            if single_value_flag:
                return values[0]
            return values

        ys += f(rs)/M

    half_light_index = np.argmin(np.abs(ys-np.max(ys)/2))
    
    """ 
    We could plot the average light distribution using:
    plt.plot(rs,ys)
    plt.show()
    It has the value in the map as a function of distance from the brightest pixel
    """
    return rs[half_light_index]

    

if __name__ == '__main__':
    """
    print(estimate_brightness(
        fits_file='Output/Maps/s13_F160w/galaxy/s13_F160w.fits',
        sigma_file='Output/Maps/s13_F160w/galaxy/s13_F160w_sigma.fits',
        mask_file='Output/Maps/s13_F160w/galaxy/s13_F160w_mask.fits',
        brightness_mask_file='Output/Maps/s13_F160w/galaxy/s13_F160w_brightness_mask.fits'
    ))
    

    print(get_pixel_in_map_from_coordinates(
        map='Output/Maps/s13_F160w/galaxy/s13_F160w.fits',
        extension=0,
        ra = 150.07597,
        dec = 2.2118269
    ))
    """
    estimate_half_light_radius(
        map='Output/Maps/s13_F160w/galaxy/s13_F160w.fits',
        extension=0,
        brightness_mask='Output/Maps/s13_F160w/galaxy/s13_F160w_brightness_mask.fits'
    )