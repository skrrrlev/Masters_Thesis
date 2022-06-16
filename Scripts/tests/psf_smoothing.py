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
from matplotlib import pyplot as plt, tight_layout
from os import makedirs, getcwd
from os.path import isdir, isfile
from shutil import rmtree, copyfile, move
from astropy.io import fits
from astropy.modeling.models import Sersic2D
from astropy.visualization import simple_norm
from scipy.signal import convolve

from MTLib.fitstools import create_psf, create_from
from MTLib.galfit import create_sersic, create_input_file, run_galfit
from MTLib import files

class TestPsfSmoothing:

    def __init__(self, test_path:str, original_psf_file:str, rel_map_size: int) -> None:
        self.test_path = test_path + '/'*(not test_path.endswith('/'))
        self.original_psf_file = original_psf_file

        with fits.open(original_psf_file) as hdul:
            data: np.ndarray = hdul[0].data
            self.psf_map = np.copy(data)
            header = hdul[0].header

        deg_per_pixel_y = float(header['CD1_1'])
        self.dy = 3600*deg_per_pixel_y
        deg_per_pixel_x = float(header['CD2_2'])
        self.dx = 3600*deg_per_pixel_x

        psf_pixel_height, psf_pixel_width = self.psf_map.shape # height is number of rows, and width is number of columns
        self.shape = (psf_pixel_height*rel_map_size, psf_pixel_width*rel_map_size)

        self.models = {}
        self.results: dict[str,dict[str,list[float]]] = {}
        self.accuracy: dict[str,dict[str,list[float]]] = {}
        self.smoothing = []

    def add_model(self, name:str, effective_radius: float, sersic_index: float, position_angle:float, axis_ratio:float, signal_to_noise_ratio:float):
        path = f'{self.test_path}{name}'
        if not isdir(path):
            makedirs(path)
        self.models[name] = {
            're':effective_radius,
            'n':sersic_index,
            'pa':position_angle,
            'ar':axis_ratio,
            'snr':signal_to_noise_ratio,
            'path':path,
            'model':f'{path}/model.fits',
            'convolved':f'{path}/convolved.fits',
            'input':f'{path}/input.fits',
            'sigma':f'{path}/sigma.fits',
            'out':{}
        }
        self.results[name] = {'re':[],'n':[],'pa':[],'ar':[],'mag':[],'δre':[],'δn':[],'δpa':[],'δar':[],'δmag':[]}
        self.accuracy[name] = {'re':[],'n':[],'pa':[],'ar':[],'mag':[],'δre':[],'δn':[],'δpa':[],'δar':[],'δmag':[]}

    def create_synthetic_maps(self):
        for name in self.models:
            
            # Create the basic model
            create_sersic(
                filename=self.models[name]["model"],
                re=self.models[name]['re'],
                n=self.models[name]['n'],
                shape=self.shape,
                dx=self.dx,
                dy=self.dy,
                **{'position_angle':self.models[name]['pa'],'axis_ratio':self.models[name]['ar']}
            )

            # Convolve the model with the psf
            with create_from(from_file=self.models[name]["model"],to_file=self.models[name]["convolved"],index=0) as hdu:
                hdu.data = convolve(hdu.data, self.psf_map, mode='same')
                peak_brightness = np.max(hdu.data)

            # Create the noise distribution
            std = peak_brightness/self.models[name]['snr']
            noise = np.random.normal(loc=0,scale=std,size=self.shape)

            # Add the noise to the convolved model
            with create_from(from_file=self.models[name]['convolved'], to_file=self.models[name]['input'], index=0) as hdu:
                hdu.data += noise

            # Create the sigma-map
            with create_from(from_file=self.models[name]['input'],to_file=self.models[name]["sigma"], index=0) as hdu:
                hdu.data = np.ones(self.shape)*std
            
    def fit_maps(self,smoothing:"list[int]"):
        self.smoothing = smoothing
        for name in self.models:
            out_path = f'{self.models[name]["path"]}/fits'
            if not isdir(out_path):
                makedirs(out_path)

            for s in smoothing:
                
                # Create the galfit input file
                self.models[name]['out'][s] = f'{out_path}/fit_s{s}.fits'
                parameters = {
                    'A': [self.models[name]['input']],
                    'B': [self.models[name]['out'][s]],
                    'C': [self.models[name]['sigma']],
                    'D': [f'{self.test_path}/psfs/psf_s{s}.fits'],
                    'J': ['26.563'],
                    'K': [str(self.dx), str(self.dy)],
                } 
                pixel_re = self.models[name]['re']/self.dx
                sersic = {
                    '0': ['sersic'],
                    '1': [f'{self.shape[1]/2+1}',f'{self.shape[0]/2+1}','0','0'],
                    '3': ['20.0','1'],
                    '4': [f'{pixel_re}','1'],
                    '5': [f'{self.models[name]["n"]}','1'],
                    '9': [f'{self.models[name]["ar"]}','1'],
                    '10': [f'{self.models[name]["pa"]}','1'],
                    'Z': ['0']
                }

                sky_bg = {
                    '0': ['sky'],
                    '1': ['0','1'],
                    '2': ['0.0000','0'],
                    '3': ['0.0000','0'],
                    'Z': ['0']
                }

                input_file = f'{self.models[name]["path"]}/galfit.input'
                create_input_file(input_file, parameters, [sersic,sky_bg])
                run_galfit(input_file, log_file=self.models[name]['out'][s].replace('.fits','.log'))

    def find_maps(self,smoothing:"list[int]"):
        self.smoothing = smoothing
        for name in self.models:
            out_path = f'{self.models[name]["path"]}/fits'
            if not isdir(out_path):
                print(f'Could not find any fits for model: "{name}".')
                continue
            for s in smoothing:
                current_file = f'{out_path}/fit_s{s}.fits'
                if isfile(current_file):
                    self.models[name]['out'][s] = f'{out_path}/fit_s{s}.fits'
                else:
                    raise ValueError(f'No such file as {current_file}.')

    def read_result_from_fit(self,name,key):
        with fits.open(self.models[name]['out'][key]) as hdul:
            header = hdul[2].header
            
            fit_string = header['1_RE']
            if '*' in fit_string:
                fit_string = fit_string.replace('*','')
                print(f'Uncertain fit of "effective radius" for {name} with smoothing {key}')
            fit_val, fit_error = findall(r'([0-9]+\.[0-9]+) \+/- ([0-9]+\.[0-9]+)', fit_string)[0]
            self.results[name]['re'].append(float(fit_val))
            self.results[name]['δre'].append(float(fit_error))
            
            fit_string = header['1_MAG']
            if '*' in fit_string:
                fit_string = fit_string.replace('*','')
                print(f'Uncertain fit of "integrated magnitude" for {name} with smoothing {key}')
            fit_val, fit_error = findall(r'([0-9]+\.[0-9]+) \+/- ([0-9]+\.[0-9]+)', fit_string)[0]
            self.results[name]['mag'].append(float(fit_val))
            self.results[name]['δmag'].append(float(fit_error))
            
            fit_string = header['1_N']
            if '*' in fit_string:
                fit_string = fit_string.replace('*','')
                print(f'Uncertain fit of "sersic index" for {name} with smoothing {key}')
            fit_val, fit_error = findall(r'([0-9]+\.[0-9]+) \+/- ([0-9]+\.[0-9]+)', fit_string)[0]
            self.results[name]['n'].append(float(fit_val))
            self.results[name]['δn'].append(float(fit_error))
            
            fit_string = header['1_AR']
            if '*' in fit_string:
                fit_string = fit_string.replace('*','')
                print(f'Uncertain fit of "effective radius" for {name} with smoothing {key}')
            fit_val, fit_error = findall(r'([0-9]+\.[0-9]+) \+/- ([0-9]+\.[0-9]+)', fit_string)[0]
            self.results[name]['ar'].append(float(fit_val))
            self.results[name]['δar'].append(float(fit_error))
            
            fit_string = header['1_PA']
            if '*' in fit_string:
                fit_string = fit_string.replace('*','')
                print(f'Uncertain fit of "effective radius" for {name} with smoothing {key}')
            fit_val, fit_error = findall(r'([0-9]+\.[0-9]+) \+/- ([0-9]+\.[0-9]+)', fit_string)[0]
            self.results[name]['pa'].append(float(fit_val))
            self.results[name]['δpa'].append(float(fit_error))
    
    def calculate_accuracy(self):
        for name in self.models:
            for result in ['re','n','mag','pa','ar']:
                fitted_values = self.results[name][result]
                fitted_errors = self.results[name][f'δ{result}']
                if result == 'mag':
                    actual_val = 20.0
                elif result == 're':
                    actual_val = self.models[name][result]/self.dx
                else:
                    actual_val = self.models[name][result]
                self.accuracy[name][result] = (np.array(fitted_values)-actual_val)/actual_val
                self.accuracy[name][f'δ{result}'] = np.array(fitted_errors)/(actual_val*np.array(fitted_values))
    
    def read_results(self):
        for name in self.models:
            for key in self.models[name]['out']:
                self.read_result_from_fit(name,key)
            
    def get_results(self):
        return self.results
            
    def get_accuracy(self):
        return self.accuracy

    @staticmethod
    def plot(smoothing:"list[int]", result:dict, percentage=False):

        number_of_models = len(result)
        number_of_cols = int(np.ceil(np.sqrt(number_of_models)))
        number_of_rows = int(np.ceil(number_of_models/number_of_cols))

        fig, axes = plt.subplots(number_of_rows,number_of_cols,sharex=True)
        fig.set_size_inches(3*number_of_cols,3*number_of_rows, forward=True)
        #plt.subplots_adjust(wspace=0)
        flat_axes = axes.flatten()
        if percentage:
            factor=100
        else:
            factor=1
        for i, name in enumerate(result):
            ax: plt.axes._subplots.AxesSubplot = flat_axes[i]
            ax.set_title(name.replace('-0','').replace('-', ' '))
            
            ax.plot(smoothing, factor*result[name]['re'], f'-k',label='effective radius',zorder=5)
            ax.scatter(smoothing,factor*result[name]['re'],color='k',s=10,zorder=10)
            ax.plot(smoothing, factor*result[name]['n'], f'--k',label='sersic index',zorder=5)
            ax.scatter(smoothing,factor*result[name]['n'],color='k',s=10,zorder=10)
            ax.legend(frameon=False,prop={'size': 7})
            ax.grid(alpha=0.25)
            if i>=number_of_cols*(number_of_rows-1):
                ax.set_xlabel('Smoothing factor')

        # Handle empty axis
        for j in range(number_of_models, len(flat_axes)):
            fig.delaxes(flat_axes[j])

        return fig, axes

    @staticmethod
    def add_to_plot(axes, smoothing, result, j, percentage=False):
        flat_axes = axes.flatten()
        if percentage:
            factor=100
        else:
            factor=1
        for i, name in enumerate(result):
            ax: plt.axes._subplots.AxesSubplot = flat_axes[i]           
            ax.plot(smoothing, factor*result[name]['re'], f'-C{j}',label='effective radius',zorder=5)
            ax.scatter(smoothing,factor*result[name]['re'],color=f'C{j}',s=10,zorder=10)
            ax.plot(smoothing, factor*result[name]['n'], f'C{j}--',label='sersic index',zorder=5)
            ax.scatter(smoothing,factor*result[name]['n'],color=f'C{j}',s=10,zorder=10)

# Constants
FIT = False
REL_MAP_SIZE = 2    # Relative size of the sersic map to the PSF map.
PSF_FWHM = 0.25     # arcseconds
RE_SMALL, RE_MEDIUM, RE_LARGE = PSF_FWHM*0.8/2, PSF_FWHM*3/2, PSF_FWHM*10/2
SNR_LOW, SNR_HIGH = 25, 100

SMOOTHING_LIST = [2,3,4,5,6,7,8]

np.random.seed(1)

# Input
psf_source_file = 'Output/Maps/HST/source20/psf/source20_F160w.fits'
psf_sigma_file = 'Output/Maps/HST/source20/psf/source20_F160w_sigma.fits'
source_region_file = 'Data/Maps/HST/source20/regions/source20_F160w.reg'

def setup():
    test_name = f'psf_smoothing_test_{files.extract_filename(psf_source_file)}'
    test_path = f'Output/tests/{test_name}/'
    if FIT:
        if isdir(test_path):
            rmtree(test_path)
        makedirs(test_path)
        copyfile(psf_source_file, f'{test_path}/source.fits')
        copyfile(psf_sigma_file, f'{test_path}/sigma.fits')
        copyfile(source_region_file, f'{test_path}/region.reg')

    return test_path

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
    copyfile(psf_file,f'{test_path}/original_psf.fits')

    if FIT:
        makedirs(f'{test_path}/psfs/')
        for s in SMOOTHING_LIST:
            create_psf(
                fits_file=f'{test_path}/source.fits',
                sigma_file=f'{test_path}/sigma.fits',
                region_file=f'{test_path}/region.reg',
                smooth_size=s
            )
            move(f'{test_path}/source_psf.fits',f'{test_path}/psfs/psf_s{s}.fits')

    accuracies = {}
    for i in range(10):
        ctrl = TestPsfSmoothing(test_path, f'{test_path}/original_psf.fits', REL_MAP_SIZE)
        re_offset = (1.2 - 0.8) * np.random.random_sample() + 0.8
        n_offset = (1.2 - 0.8) * np.random.random_sample() + 0.8

        ctrl.add_model(f'unresolved-low-snr-{i}', RE_SMALL*re_offset, 4*n_offset, 30, 0.8, SNR_LOW)
        ctrl.add_model(f'resolved-low-snr-{i}', RE_MEDIUM*re_offset, 2*n_offset, 30, 0.8, SNR_LOW)
        ctrl.add_model(f'extended-low-snr-{i}', RE_LARGE*re_offset, 2*n_offset, 30, 0.8, SNR_LOW)

        ctrl.add_model(f'unresolved-high-snr-{i}', RE_SMALL*re_offset, 4*n_offset, 30, 0.8, SNR_HIGH)
        ctrl.add_model(f'resolved-high-snr-{i}', RE_MEDIUM*re_offset, 2*n_offset, 30, 0.8, SNR_HIGH)
        ctrl.add_model(f'extended-high-snr-{i}', RE_LARGE*re_offset, 2*n_offset, 30, 0.8, SNR_HIGH)

        ctrl.create_synthetic_maps()
        if FIT:
            ctrl.fit_maps(SMOOTHING_LIST)
        else:
            ctrl.find_maps(SMOOTHING_LIST)

        ctrl.read_results()
        ctrl.calculate_accuracy()
        accuracy = ctrl.get_accuracy()
        accuracies[i] = accuracy

        if not i:
            fig, axes = TestPsfSmoothing.plot(SMOOTHING_LIST, accuracy, percentage=True)
            for ax in axes[:,0]:
                ax.set_ylabel('Percentage Residuals')
        else:
            TestPsfSmoothing.add_to_plot(axes,SMOOTHING_LIST,accuracy,i-1,percentage=True)

    plt.savefig('figures/psfs/combined_smoothing_test_residuals.png',bbox_inches='tight')

    '''
    for name in ctrl.models:
        print(f'Model: "{name}"')
        print(f're = {ctrl.models[name]["re"]/ctrl.dx:.3f}, n = {ctrl.models[name]["n"]:.3f}, mag = 20.000, pa = {ctrl.models[name]["pa"]:.3f}, ar = {ctrl.models[name]["ar"]:.3f}')
        for key in ['re','n','mag','pa','ar']:
            print(f'\t{key} = ')
            print(f'\t\tAccuracy:  {", ".join([f"{val:.3f}" for val in accuracy[name][key]])}')

    
    fig, axes = ctrl.plot(accuracy, percentage=True)
    for ax in axes[:,0]:
        ax.set_ylabel('Percentage Residuals')
    plt.savefig('figures/psfs/smoothing_test_residuals.png',bbox_inches='tight')
    '''
    for key in accuracies:
        if not key:
            accuracy: dict = accuracies[key]
            names = list(accuracy.keys())
            for name in names:   
                accuracy[name]['re'] /= 10
                accuracy[name]['n'] /= 10
            continue
        
        other:dict = accuracies[key]
        other_names = list(other.keys())
        for name,other_name in zip(names,other_names):
            accuracy[name]['re'] += other[other_name]['re']/10
            accuracy[name]['n'] += other[other_name]['n']/10
    
    fig, axes = TestPsfSmoothing.plot(SMOOTHING_LIST, accuracy, percentage=True)
    for ax in axes[:,0]:
        ax.set_ylabel('Percentage Residuals')
    plt.savefig('figures/psfs/average_smoothing_test_residuals.png',bbox_inches='tight')

if __name__ == '__main__':
    main()