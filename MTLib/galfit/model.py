
"""
# IMAGE PARAMETERS
A) none                         # Input data image (FITS file)
B) {filename}                   # Output data image block
C) none                         # Sigma image name (made from data if blank or "none")
D) none                         # Input PSF image and (optional) diffusion kernel
E) 1                            # PSF fine sampling factor relative to data
F) none                         # Bad pixel mask (FITS image or ASCII coord list)
G) none                         # File with parameter constraints (ASCII file)
H) 0 {shape[1]} 0 {shape[0]}    # Image region to fit (xmin xmax ymin ymax)
I) {shape[1]} {shape[0]}        # Size of the convolution box (x y)
J) {mag_zp}                     # Magnitude photometric zeropoint
K) {dx} {dy}                    # Plate scale (dx dy) [arcsec per pixel]
O) both                         # Display type (regular, curses, both)
P) 1                            # Options: 0=normal run; 1,2=make model/imgblock & quit

# INITIAL FITTING PARAMETERS #
# ------------------------------------------------------------------------------
# par) par value(s) fit toggle(s) # parameter description
# ------------------------------------------------------------------------------

0) sersic                   # object type
1) {pos_x} {pos_y} 0 0      # position x, y
3) {mag} 0                  # integrated magnitude
4) {re} 0                   # R_e (half-light radius) [pix]
5) {n} 0                    # Sersic index n (de Vaucouleurs n=4)
9) {ar} 0                   # axis ratio (b/a)
10) {pa} 0                  # position angle (PA) [deg: Up=0, Left=90]
Z) 0                        # output option (0 = resid., 1 = Don’t subtract)
"""

import numpy as np
from astropy.io import fits
from os.path import isdir
from os import makedirs
from shutil import rmtree
from subprocess import run
from ..files import extract_path

TEMP_DIR = 'tmp/'

def _create_input_file(filename, shape, mag_zp, dx, dy, pos_x, pos_y, mag, re, n, ar, pa):
    input_file_as_string = f"""# IMAGE PARAMETERS
A) none                         # Input data image (FITS file)
B) {filename}                   # Output data image block
C) none                         # Sigma image name (made from data if blank or "none")
D) none                         # Input PSF image and (optional) diffusion kernel
E) 1                            # PSF fine sampling factor relative to data
F) none                         # Bad pixel mask (FITS image or ASCII coord list)
G) none                         # File with parameter constraints (ASCII file)
H) 0 {shape[1]} 0 {shape[0]}    # Image region to fit (xmin xmax ymin ymax)
I) {shape[1]} {shape[0]}        # Size of the convolution box (x y)
J) {mag_zp}                     # Magnitude photometric zeropoint
K) {dx} {dy}                    # Plate scale (dx dy) [arcsec per pixel]
O) both                         # Display type (regular, curses, both)
P) 1                            # Options: 0=normal run; 1,2=make model/imgblock & quit

# INITIAL FITTING PARAMETERS #
# ------------------------------------------------------------------------------
# par) par value(s) fit toggle(s) # parameter description
# ------------------------------------------------------------------------------

0) sersic                   # object type
1) {pos_x} {pos_y} 0 0      # position x, y
3) {mag} 0                  # integrated magnitude
4) {re} 0                   # R_e (half-light radius) [pix]
5) {n} 0                    # Sersic index n (de Vaucouleurs n=4)
9) {ar} 0                   # axis ratio (b/a)
10) {pa} 0                  # position angle (PA) [deg: Up=0, Left=90]
Z) 0                        # output option (0 = resid., 1 = Don’t subtract)
"""
    file_path = f'{TEMP_DIR}make_sersic.input'
    with open(file_path,'w') as f:
        f.write(input_file_as_string)

    return file_path

def create_sersic(filename: str, re:float, n:float, shape:"tuple[float]", dx:float, dy:float, **kwargs):
    '''
    Using galfit, create a fits file with a 2D Sersic profile.
    
    <filename>:     The filename (and path).
    <re>:           the effective radius in arcsec
    <shape>:        Number of pixels in (height, width)
    <dx>:           arcsec per pixel along width of map
    <dy>:           arcsec per pixel along height of map

    **kwargs:
    <pos_x>:                    x-position in number of pixels of centre of synthetic source. Defaults to width/2
    <pos_y>:                    y-position in number of pixels of centre of synthetic source. Defaults to height/2
    <position_angle>:           position angle of source in degrees. Defaults to 0
    <axis_ratio>:               axis ratio (b/a). Defaults to 1.
    <magnitude_zero_point>:     zero point of magnitude scale. Defualts to 26.563.
    <integrated_magnitude>:     integrated magnitude of synthetic source. Defaults to 20.0.
    '''

    '''create temp directory and output directory'''
    if isdir(TEMP_DIR):
        rmtree(TEMP_DIR)
    makedirs(TEMP_DIR)

    path_of_model_file = extract_path(filename)
    if path_of_model_file and not isdir(path_of_model_file):
        makedirs(path_of_model_file)
    
    '''Load parameters from kwargs'''

    map_pixel_height, map_pixel_width = shape

    if 'pos_x' in kwargs:
        pos_x = kwargs['pos_x']
    else:
        pos_x = map_pixel_width/2
    if 'pos_y' in kwargs:
        pos_y = kwargs['pos_y']
    else:
        pos_y = map_pixel_height/2
    if 'position_angle' in kwargs:
        pa = kwargs['position_angle']
    else:
        pa = 0
    if 'axis_ratio' in kwargs:
        ar = kwargs['axis_ratio']
    else:
        ar = 1
    if 'magnitude_zero_point' in kwargs:
        mag_zp = kwargs['magnitude_zero_point']
    else:
        mag_zp = 26.563
    if 'integrated_magnitude' in kwargs:
        mag = kwargs['integrated_magnitude']
    else:
        mag = 20.0

    # convert re to pixels
    re = abs(re/dx)
    
    '''Run galfit'''
    galfit_input_file = _create_input_file(filename,shape,mag_zp,dx,dy,pos_x,pos_y,mag,re,n,ar,pa)
    run(['galfit',galfit_input_file])

    '''Clean up'''
    rmtree(TEMP_DIR)



if __name__ == '__main__':

    kwargs = {
        'position_angle':30,
        'axis_ratio':0.8,
    }


    create_sersic(
        filename='test/test.fits',
        re=0.16,
        n=2,
        shape=(100,100),
        dx=0.1,
        dy=-0.1,
        **kwargs)