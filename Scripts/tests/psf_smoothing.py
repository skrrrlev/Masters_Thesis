'''
Author: Ditlev Frickmann

When I create PSFs, I select a pointsource with an appropriate "psf pattern" (e.g. aery rings, cross-artifcat) around it, to represent the psf.
To avoid the noise in the original mosaic map, I mask all pixels with a SNR > 3. I smooth the mask and apply it. 
In this script, I'd like to test how the smoothing affects Galfit in recovering model parameters. 

Idea:
- Create 2D Sersic profile
- Convolve with PSF and add Gaussian noise
- Use galfit to fit the image supplying PSFs with different smoothing parameters
- Estimate the accuracy when the smoothing parameters are off of the actual value.
- Test with two different noise realisations. (One with low S/N and one with high S/N)

'''
from sys import argv
from click import style
import numpy as np
from re import findall
from matplotlib import pyplot as plt
from os import makedirs, getcwd
from os.path import isdir, isfile
from shutil import rmtree, copyfile
from astropy.io import fits
from astropy.modeling.models import Sersic2D
from astropy.visualization import simple_norm
from scipy.signal import convolve

from MTLib.fitstools import create_psf, create_from
from MTLib.galfit import create_sersic, create_input_file, run_galfit
from MTLib import files

# Constants
REL_MAP_SIZE = 2    # Relative size of the sersic map to the PSF map.

PSF_FWHM = 0.16     # arcseconds
RE_LIST = [PSF_FWHM*0.8/2,PSF_FWHM*2/2,PSF_FWHM*3/2]
RE_STRING_LIST = ['unresolved','slightly_resolved','resolved']

SERSIC_INDEX_LIST = [4,2,2] 
MODEL_KWARGS = {
        'position_angle':30,
        'axis_ratio':0.8,
    }

NOISE_SNR_LIST = [10,50] # brightest pixel of synthetic source divided by standard deviation of added noise.
NOISE_STRING_LIST = ['low_SNR','high_SNR']

SMOOTHING_LIST = [2,3,4,5,6,7,8]

COLOURS = ['#f2bC94','#30110d','#722620']
SHAPES = ['--','-.']

np.random.seed(1)

# Input
psf_source_file = 'Output/Maps/HST/source20/psf/source20_F160w.fits'
psf_sigma_file = 'Output/Maps/HST/source20/psf/source20_F160w_sigma.fits'
source_region_file = 'Data/Maps/HST/source20/regions/source20_F160w.reg'

def setup():
    test_name = f'psf_smoothing_test_{files.extract_filename(psf_source_file)}'
    test_path = f'Output/tests/{test_name}/'
    if isdir(test_path):
        rmtree(test_path)
    
    makedirs(test_path)
    copyfile(psf_source_file, f'{test_path}/source.fits')
    copyfile(psf_sigma_file, f'{test_path}/sigma.fits')
    copyfile(source_region_file, f'{test_path}/region.reg')

    makedirs(f'{test_path}/models/')
    makedirs(f'{test_path}/convolved/')
    makedirs(f'{test_path}/noisy/')
    makedirs(f'{test_path}tmp/')

    return test_path

def plot_map(map: np.ndarray):
    fig = plt.figure(figsize=(5,5), dpi=150)

    norm = simple_norm(map, 'sqrt', percent=99.9)
    plt.imshow(map, norm=norm, origin='lower', cmap='gist_heat')

    return fig

def get_sersic_map(width,height,r_eff):
    sersic = Sersic2D(
        amplitude = SOURCE_AMPLITUDE,   # surface brightness at <r_eff>
        r_eff = r_eff,                  # effective radius (half-light radius)
        n = 2,                          # sersic index
        x_0 = width/2,
        y_0 = height/2,
        ellip = 0.3,                    # ellipticity
        theta = np.pi/3,                # rotation angle in radians
    )

    x,y = np.meshgrid(np.arange(width), np.arange(height))

    return sersic(x,y)



def main():
    
    test_path = setup()

    '''Create the PSF from the source image'''
    create_psf(
        fits_file=f'{test_path}/source.fits',
        sigma_file=f'{test_path}/sigma.fits',
        region_file=f'{test_path}/region.reg',
        smooth_size=3
    )
    psf_file = f'{test_path}/source_psf.fits'
    copyfile(psf_file,f'{test_path}/original_psf.fits' )

    '''Load in the psf_file and set some constants from it.'''
    with fits.open(psf_file) as hdul:
        data: np.ndarray = hdul[0].data
        psf_map = np.copy(data)
        header = hdul[0].header

    deg_per_pixel_y = float(header['CD1_1'])
    dy = 3600*deg_per_pixel_y
    deg_per_pixel_x = float(header['CD2_2'])
    dx = 3600*deg_per_pixel_x

    psf_pixel_height, psf_pixel_width = psf_map.shape # height is number of rows, and width is number of columns
    shape = (psf_pixel_height*REL_MAP_SIZE, psf_pixel_width*REL_MAP_SIZE)

    '''Create synthetic maps'''

    for re_arcsec, sersic_index, name in zip(RE_LIST,SERSIC_INDEX_LIST, RE_STRING_LIST):
        create_sersic(
            filename=f'{test_path}/models/{name}.fits',
            re=re_arcsec,
            n=sersic_index,
            shape=shape,
            dx=dx,
            dy=dy,
            **MODEL_KWARGS)

        with create_from(f'{test_path}/models/{name}.fits',f'{test_path}/convolved/{name}.fits',index=0) as hdu:
            hdu.data = convolve(hdu.data, psf_map, mode='same')
            peak_brightness = np.max(hdu.data)

        for snr, snr_name in zip(NOISE_SNR_LIST, NOISE_STRING_LIST):
            std = peak_brightness/snr
            noise = np.random.normal(loc=0,scale=std,size=shape)
            with create_from(f'{test_path}/convolved/{name}.fits',f'{test_path}/noisy/{name}_{snr_name}.fits', index=0) as hdu:
                hdu.data += noise
            with create_from(f'{test_path}/noisy/{name}_{snr_name}.fits',f'{test_path}/noisy/{name}_{snr_name}_sigma.fits', index=0) as hdu:
                hdu.data = np.ones(shape)*std

    total_results = []
    total_errors = []
    for re, sersic_index, name in zip(RE_LIST,SERSIC_INDEX_LIST,RE_STRING_LIST):
        for snr_name in NOISE_STRING_LIST:
            results = {'Effective radius':[],'Integrated magnitude':[],'Position angle':[],'Axis ratio':[], 'Sersic index':[]}
            errors = {'Effective radius':[],'Integrated magnitude':[],'Position angle':[],'Axis ratio':[], 'Sersic index':[]}
            for smoothing_factor in SMOOTHING_LIST:
                create_psf(
                    fits_file=f'{test_path}/source.fits',
                    sigma_file=f'{test_path}/sigma.fits',
                    region_file=f'{test_path}/region.reg',
                    smooth_size=smoothing_factor
                )

                fits_file = f'{test_path}/noisy/{name}_{snr_name}.fits'
                sigma_file = f'{test_path}/noisy/{name}_{snr_name}_sigma.fits'
                out_file = f'{test_path}tmp/out.fits'
                parameters = {
                    'A': [fits_file],
                    'B': [out_file],
                    'C': [sigma_file],
                    'D': [psf_file],
                    'J': ['26.563'],
                    'K': [str(dx), str(dy)],
                } 

                sersic = {
                    '0': ['sersic'],
                    '1': [f'{shape[1]/2}',f'{shape[0]/2}','0','0'],
                    '3': ['20.0','1'],
                    '4': [f'{re/dx}','1'],
                    '5': [f'{sersic_index}','1'],
                    '9': [str(MODEL_KWARGS['axis_ratio']),'1'],
                    '10': [str(MODEL_KWARGS['position_angle']),'1'],
                    'Z': ['0']
                }

                sky_bg = {
                    '0': ['sky'],
                    '1': ['0','1'],
                    '2': ['0.0000',0],
                    '3': ['0.0000',0],
                    'Z': ['0']
                }

                create_input_file('galfit.INPUT', parameters, [sersic,sky_bg])
                run_galfit('galfit.INPUT', log_file=f'{test_path}tmp/out.log')

                with fits.open(out_file) as hdul:
                    header = hdul[2].header
                    
                    fit_string = header['1_RE']
                    fit_val, fit_error = findall(r'([0-9]+\.[0-9]+) \+/- ([0-9]+\.[0-9]+)', fit_string)[0]
                    results['Effective radius'] += [(float(fit_val)-(re/dx))/(re/dx)]
                    errors['Effective radius'] += [abs(float(fit_error)/((re/dx)*float(fit_val)))]
                    
                    fit_string = header['1_MAG']
                    fit_val, fit_error = findall(r'([0-9]+\.[0-9]+) \+/- ([0-9]+\.[0-9]+)', fit_string)[0]
                    results['Integrated magnitude'] += [(float(fit_val)-(20.0))/(20.0)]
                    errors['Integrated magnitude'] += [abs(float(fit_error)/(20.0*float(fit_val)))]
                    
                    fit_string = header['1_AR']
                    fit_val, fit_error = findall(r'([0-9]+\.[0-9]+) \+/- ([0-9]+\.[0-9]+)', fit_string)[0]
                    results['Position angle'] += [(float(fit_val)-(MODEL_KWARGS['axis_ratio']))/(MODEL_KWARGS['axis_ratio'])]
                    errors['Position angle'] += [abs(float(fit_error)/(MODEL_KWARGS['axis_ratio']*float(fit_val)))]
                    
                    fit_string = header['1_PA']
                    fit_val, fit_error = findall(r'([0-9]+\.[0-9]+) \+/- ([0-9]+\.[0-9]+)', fit_string)[0]
                    results['Axis ratio'] += [(float(fit_val)-(MODEL_KWARGS['axis_ratio']))/(MODEL_KWARGS['axis_ratio'])]
                    errors['Axis ratio'] += [abs(float(fit_error)/(MODEL_KWARGS['axis_ratio']*float(fit_val)))]

                    fit_string:str = header['1_N']
                    if '*' in fit_string:
                        fit_string = fit_string.replace('*','')
                        print(f'Error in "{name} {snr_name}".')
                    fit_val, fit_error = findall(r'([0-9]+\.[0-9]+) \+/- ([0-9]+\.[0-9]+)', fit_string)[0]
                    results['Sersic index'] += [(float(fit_val)-(sersic_index))/(sersic_index)]
                    errors['Sersic index'] += [abs(float(fit_error)/(sersic_index*float(fit_val)))]

            total_results.append(results)
            total_errors.append(errors)
    
    print(total_results)
    print()
    print(total_errors)

    for key in ['Effective radius','Integrated magnitude','Position angle','Axis ratio', 'Sersic index']:
        plt.figure(figsize=(5,5),dpi=100)
        plt.title(key)
        i = 0
        for name,colour in zip(RE_STRING_LIST,COLOURS):
            for snr_name,shape in zip(NOISE_STRING_LIST,SHAPES):
                plt.plot(SMOOTHING_LIST,total_results[i][key],shape,color=colour,label=f'{name} {snr_name}')
                i += 1
        plt.legend(frameon=False)
        plt.savefig(f'par_{key.replace(" ", "_")}.png')

if __name__ == '__main__':
    main()