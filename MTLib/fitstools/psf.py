import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, FK5, ICRS
from astropy.visualization import simple_norm
import matplotlib.pyplot as plt
from re import findall
from scipy.ndimage import median_filter
from os.path import isdir
from os import makedirs

from .. import files, constants

def _mask_psf(data: np.ndarray, sigma: np.ndarray, smooth_size:int=3):
    '''Create and apply the sigma mask for the PSF'''
    mask = data/sigma < 3
    smoothed_mask = median_filter(mask,size=smooth_size)
    data[smoothed_mask] = 0

    return data

def _parse_region_file(region_file:str):
    with open(region_file,'r') as rf:
        content = rf.readlines()
    
    regions: dict[str,dict[str,float]] = {}
    for line in content:
        matches = findall(r'(\d+(?:\.\d+)?(?:e[+-]?\d+)?)\S', line)
        if len(matches) >= 5:
            try:
                region_name:str = findall(r'text=\{(\S+)\}', line)[0].lower()
            except IndexError: # A name was not defined for the region
                continue
            xcentre:float = float(matches[0])-1
            ycentre:float = float(matches[1])-1
            xwidth:float = float(matches[2])
            ywidth:float = float(matches[3])

            xmax = xcentre+xwidth/2
            xmin = xcentre-xwidth/2

            ymin = ycentre-ywidth/2
            ymax = ycentre+ywidth/2
            regions[region_name] = {'bottom-left':{'x':xmin,'y':ymin},'top-right':{'x':xmax,'y':ymax}}
            # print(f'{region_name}: {xmin}, {xmax}, {ymin}, {ymax}')
    return regions

def _apply_region_mask(data:np.ndarray, region_file:str):
    '''If a region with the name "out" is overlapping the psf region, set all the values to zero in the psf for that region'''
    
    regions = _parse_region_file(region_file)
    if 'out' in regions:
        out = regions.pop('out')
    else:
        return data
    
    for region_name in regions:
        region = regions[region_name]
        if out['bottom-left']['x']>=region['top-right']['x'] or region['bottom-left']['x']>=out['top-right']['x']or out['top-right']['y']<=region['bottom-left']['y'] or region['top-right']['y']<=out['bottom-left']['y']:
            # out does not overlap with the region
            continue

        ylim = region['top-right']['y']-region['bottom-left']['y']
        ymin = int(max(0,np.floor((out['bottom-left']['y']-region['bottom-left']['y']))))
        ymax = int(min(ylim,np.ceil(out['top-right']['y']-region['bottom-left']['y']-1)))
        xlim = region['top-right']['x']-region['bottom-left']['x']
        xmin = int(max(0,np.floor(out['bottom-left']['x']-region['bottom-left']['x'])))
        xmax = int(min(xlim,np.ceil(out['top-right']['x']-region['bottom-left']['x'])))

        data[ymin:ymax,xmin:xmax] = 0

    return data

def _normalize_psf(data: np.ndarray):
    '''The sum of all pixels is 1.'''
    return data/np.sum(data)

class CenteringError(Exception):
    '''Raised when the centering fails'''

def _centre_psf(data: np.ndarray):
    '''Concatenate blank rows and columns to the data array, until the PSF is centered correctly.'''

    row_flag, col_flag = False, False
    i = 0
    while True:

        '''Extract the centre of the map as the place with the brightest pixel.'''
        centre = np.array(np.where(data==np.amax(data))).flatten()
        row_centre_index = centre[0]
        col_centre_index = centre[1]
        N,M = data.shape

        '''
        In the following code, rows and columns are concatenated in the correct places, until the specific condition is achieved. 
            odd amount of pixels: centre should be at NPIX/2+0.5
            even amount of pixels: centre should be at NPIX/2+1
        '''

        '''Treat the rows'''
        if N%2: # odd: centre should be at NPIX/2+0.5
            desired_row_index = N/2+0.5 -1
        else: # even: centre should be at NPIX/2+1
            desired_row_index = N/2+1 -1  
        if row_centre_index<desired_row_index: # actual centre lies above desired centre
            data = np.concatenate((np.zeros((1,M)),data),axis=0) # add empty line on top of map
            N += 1
        elif row_centre_index>desired_row_index: # actual centre lies below desired centre
            data = np.concatenate((data,np.zeros((1,M))),axis=0) # add empty line on bottom of map
            N += 1
        else:
            row_flag = True

        '''Treat the columns'''
        if M%2: # odd: centre should be at NPIX/2+0.5
            desired_col_index = M/2+0.5 -1
        else: # even: centre should be at NPIX/2+1
            desired_col_index = M/2+1 -1
        if col_centre_index>desired_col_index: # actual centre lies right of desired centre
            data = np.concatenate((data,np.zeros((N,1))),axis=1) # add empty column on rhs of map
        elif col_centre_index<desired_col_index: # actual centre is left of desired centre
            data = np.concatenate((np.zeros((N,1)),data),axis=1) # add empty column on lhs of map
        else:
            col_flag = True
        
        if row_flag and col_flag:
            #print(f'Shape: {data.shape}')
            #print(f'Centre: {np.array(np.where(data==np.amax(data))).flatten()+1}')
            break
        else:
            i += 1
            if i > 15:
                raise CenteringError(f'Could not centre PSF in 15 iterations of the centering algorithm. Current shape: {data.shape}, Centre: {np.array(np.where(data==np.amax(data))).flatten()+1}') 

    return data

def _plot(data: np.ndarray, filename:str, adjust:str):
    '''2d image plot with profiles. Saved using the <filename>. <adjust> must be either "left" or "right".'''

    h, w = data.shape
    x = np.arange(w)
    y = np.arange(h)
    X,Y = np.meshgrid(x,y)
    data_ratio = h / float(w)

    # Get the profile along the centre of the map
    if w%2:
        W = w/2+0.5 -1
    else:
        W = w/2+1 -1

    if h%2:
        H = h/2+0.5 -1
    else:
        H = h/2+1 -1

    ''' define axes lenght in inches '''

    width_ax0 = 8.
    width_ax1 = 2.
    height_ax2 = 2.

    height_ax0 = width_ax0 * data_ratio

    ''' define margins size in inches '''

    left_margin  = 0.65
    right_margin = 0.2
    bottom_margin = 0.5
    top_margin = 0.25
    inter_margin = 0.5

    ''' calculate total figure size in inches '''

    fwidth = left_margin + right_margin + inter_margin + width_ax0 + width_ax1
    fheight = bottom_margin + top_margin + inter_margin + height_ax0 + height_ax2

    fig = plt.figure(figsize=(fwidth, fheight))
    fig.patch.set_facecolor('white')

    '''create axes'''
    if adjust == 'right':
        ax0 = fig.add_axes([left_margin / fwidth,
                            (bottom_margin + inter_margin + height_ax2) / fheight,
                            width_ax0 / fwidth, height_ax0 / fheight])

        ax1 = fig.add_axes([(left_margin + width_ax0 + inter_margin) / fwidth,
                            (bottom_margin + inter_margin + height_ax2) / fheight,
                            width_ax1 / fwidth, height_ax0 / fheight])

        ax2 = fig.add_axes([left_margin / fwidth, bottom_margin / fheight,
                            width_ax0 / fwidth, height_ax2 / fheight])     
    elif adjust == 'left':
        ax0 = fig.add_axes([(left_margin + width_ax1 + inter_margin) / fwidth,
                            (bottom_margin + inter_margin + height_ax2) / fheight,
                            width_ax0 / fwidth, height_ax0 / fheight])

        ax1 = fig.add_axes([left_margin / fwidth,
                            (bottom_margin + inter_margin + height_ax2) / fheight,
                            width_ax1 / fwidth, height_ax0 / fheight])

        ax2 = fig.add_axes([(left_margin + width_ax1 + inter_margin) / fwidth, bottom_margin / fheight,
                            width_ax0 / fwidth, height_ax2 / fheight])
    else:
        raise ValueError('<adjust> must be "left" or "right".')
    
    '''plot data'''

    bounds = [x.min(),x.max(),y.min(),y.max()]
    norm = simple_norm(data, 'sqrt', percent=99.)
    ax0.imshow(data, norm=norm, origin='lower', cmap='gist_heat',extent = bounds)
    ax0.xaxis.set_ticklabels([])
    ax0.yaxis.set_ticklabels([])
    ax0.grid(alpha=0.25)

    ax1.plot(data[:,int(W)],Y[:,int(W)],'k.',data[:,int(W)],Y[:,int(W)],'k-')
    ax1.grid(alpha=0.25)
    if adjust == 'left':
        ax1.invert_xaxis()
    ax1.set_xlabel('Unit of map')
    ax1.set_ylim((y.min(),y.max()))
    ax1.yaxis.set_ticklabels([])
    
    ax2.plot(X[int(H),:], data[int(H),:], 'k.', X[int(H),:], data[int(H),:],'k-')
    ax2.grid(alpha=0.25)
    ax2.set_xlim((x.min(),x.max()))
    ax2.xaxis.set_ticklabels([])
    ax2.set_ylabel('Unit of map')
    
    path = files.extract_path(filename)
    if not isdir(path):
        makedirs(path)
    fig.savefig(filename,bbox_inches = 'tight')
    plt.close()

def create_psf(fits_file:str, sigma_file:str, region_file:str='', smooth_size=3, output_file:str = ''):
    ''''''
    with fits.open(fits_file) as hdul:
        data: np.ndarray = hdul[0].data
        psf = np.copy(data)
    with fits.open(sigma_file) as hdul:
        sigma: np.ndarray = hdul[1].data
    
    psf = _mask_psf(psf, sigma, smooth_size)
    if region_file:
        psf = _apply_region_mask(psf, region_file)
    psf = _normalize_psf(psf)
    psf = _centre_psf(psf)
    #plt.show()

    '''Save the psf.'''
    # Get the weight image in the fits file
    with fits.open(fits_file) as hdul:
        hdu = hdul[0]
        hdu.header['EXTVER'] = 'PSF'
        hdu.data = psf

        # Now we save the psf
        if not output_file:
            current_filename = files.extract_filename(fits_file)
            output_folder = files.extract_path(fits_file).replace('Data','Output')
            output_name = current_filename + '_psf'
            output_file = output_folder + output_name + '.fits'
        else:
            output_name = files.extract_filename(output_file)
        print(f'Saving psf: {output_file}')
        hdu.writeto(output_file, overwrite=True)

    _plot(psf,f'{constants.PATH.FIGURES.value}psfs/{output_name}.png',adjust="right")
    _plot(data,f'{constants.PATH.FIGURES.value}psfs/{output_name}_before_psf.png',adjust="left")

if __name__ == '__main__':
    create_psf(
        'Output/Maps/HST/source20/psf/source20_F160w.fits',
        'Output/Maps/HST/source20/psf/source20_F160w_sigma.fits',
        'Data/Maps/HST/source20/regions/source20_F160w.reg'
    )