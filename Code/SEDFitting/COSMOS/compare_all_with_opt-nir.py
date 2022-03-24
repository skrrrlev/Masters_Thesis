# imports
import sys
from os.path import isfile
from astropy.io import fits
from typing import NoReturn, Any
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
import numpy as np

# MTLib
from MTLib import cosmos

# definitions
files: "list[str]" = []

# Other methods
def usage() -> NoReturn:
    usage = '''
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    usage: python __file__.py output1 output2
    \tRun from the root directory of this project.\n
    Parameters:
               output: the output FITS file of a Stardust run.
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    '''
    print(usage)
    exit()

def setup():
    number_of_input = len(sys.argv) - 1
    if not number_of_input == 2:
        usage()
    
    for input_index in range(1, number_of_input+1):
        file = sys.argv[input_index]

        if file.endswith('.fits') and isfile(file):
            files.append(file)
        else:
            raise ValueError(f'Input {input_index} ({file}) is not a valid FITS file.')
    



# main script

def main():
    setup()
    # Load data from each file
    print('Loading data ...')
    data: "dict[str,dict[int,Any]]" = {}
    ids: "list[int]" = []
    for file in files:
        print(f'\tLoading from {file}')
        file_name = file.replace('\\','/').split('/')[-1]
        data[file_name] = {}
        with fits.open(file,'readonly') as f:
            HDU: fits.hdu.table.BinTableHDU = f[1]
            fdata: fits.fitsrec.FITS_rec = HDU.data
        for i, idx in enumerate(fdata.field('id')):
            data[file_name][idx] = fdata.field('mstar')[i]
            if idx not in ids:
                ids.append(idx)
    file_names = [file.strip('.fits') for file in list(data.keys())]
    print(data)
    print('... Done!')
    
    # Load catalogue value
    print('Loading catalogue masses... ')
    cat = cosmos.Catalogue()
    cat_mass = {idx:float(cat.get_property([idx], 'ez_mass')[0]) for idx in ids}
    print(cat_mass)
    print('... Done!')

    print('Creating figure ...')
    plt.figure(num=1, figsize=(7,5), dpi=150)
    rot = [0,45]
    y = np.arange(len(ids))+1
    for i,file in enumerate(data):
        print(f'Adding from {file}')
        x = np.array([data[file][idx]/cat_mass[idx] for idx in ids])
        t = MarkerStyle(marker='x')
        t._transform = t.get_transform().rotate_deg(rot[i])
        plt.scatter(x,y,label=file, marker=t, c='k',alpha=0.8)
            
    plt.legend(frameon=False)
    plt.grid(alpha=0.2)
    plt.xscale('log')
    ylim = [0,y[-1]+1]
    plt.ylim(ylim)
    plt.plot([1,1],ylim,'k',alpha=0.8)
    plt.yticks(y,[f'ID{idx}' for idx in ids])
    plt.xlabel(r'Fitted stellar mass / COSMOS stellar mass')
    plt.savefig(f'compare_{file_names[0]}_{file_names[1]}.png')
    print('... Done!')


if __name__ == '__main__':
    main()
    exit()
