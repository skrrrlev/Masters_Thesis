"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

This script plots the SED of a range of temperatures given by the <temperatures> array.
Beneath, Wiens displacement law is plotted for the same wavelengths.

"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
from astropy import constants as c
from astropy import units as u


planck =  lambda ν,T: (((2*c.h*ν**3)/c.c**2) * 1/(np.exp((c.h*ν)/(c.k_B*T))-1)).to(u.Jy)


def main():
    temperatures = np.arange(1,10,2,dtype=int)*1000*u.K
    N = 1000
    wavelengths = np.logspace(-2,2,N)*u.um
    frequencies = (c.c / wavelengths).to(u.GHz)
    
    fig = plt.figure()
    gs1 = GridSpec(2, 1, height_ratios=[5, 1], hspace=0.05)
    ax1 = fig.add_subplot(gs1[0])
    ax2 = fig.add_subplot(gs1[1])
    print(type(ax1))
    colours = ['red','orange','yellow','green','blue']

    for t,colour in zip(temperatures,colours):
        ax1.plot(wavelengths.value,frequencies.value*planck(frequencies,t).to(u.uJy).value,label=f'{int(t.value)}K',c=colour)
    
    ax1.set_ylabel(r'Reduced brightness $\nu I_\nu$ [GHz Jy]')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylim((1E9,1E32))
    ax1.set_xlim((1E-2,1E2))
    ax1.grid(alpha=0.25)
    ax1.legend(frameon=False)
    ax1.set_xticklabels([])

    T = (c.b_wien/wavelengths).to(u.K)
    ax2.plot(wavelengths.value,T.value, 'k')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.grid(alpha=0.25)
    ax2.set_xlabel(r'Wavelength [$\mu m$]')
    ax2.set_ylabel(r'Temperature [K]')
    ax2.set_xlim((1E-2,1E2))
    ax2.set_ylim((5E2,5E4))

    plt.savefig('bb.png')

if __name__ == '__main__':
    main()