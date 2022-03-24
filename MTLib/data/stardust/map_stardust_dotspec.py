# imports
from dataclasses import dataclass, field
from astropy import units as U

# relative imports
from .dataclasses import Filter

def get_cosmos_dotspec_map():
    '''Returns a list of filters that are in the correct order of the .spec files supplied for the COSMOS survey.'''
    stardust_dotspec_map = [
        Filter(
            instrument = 'MegaCam',
            telescope = 'CFHT',
            band = 'u*',
            spec_wavelength = 3896*U.AA,
            spec_width = 597.6*U.AA,
            stardust_code = 352
        ),
        Filter(
            instrument = 'MegaCam',
            telescope = 'CFHT',
            band = 'u',
            spec_wavelength = 3739*U.AA,
            spec_width = 518.1*U.AA,
            stardust_code = 353
        ),
        Filter(
            instrument = 'HSC',
            telescope = 'Subaru',
            band = 'g',
            spec_wavelength = 4724*U.AA,
            spec_width = 1387*U.AA,
            stardust_code = 314
        ),
        Filter(
            instrument = 'HSC',
            telescope = 'Subaru',
            band = 'r',
            spec_wavelength = 6110*U.AA,
            spec_width = 1546*U.AA,
            stardust_code = 315
        ),
        Filter(
            instrument = 'HSC',
            telescope = 'Subaru',
            band = 'i',
            spec_wavelength = 7612*U.AA,
            spec_width = 1469*U.AA,
            stardust_code = 316
        ),
        Filter(
            instrument = 'HSC',
            telescope = 'Subaru',
            band = 'z',
            spec_wavelength = 8901*U.AA,
            spec_width = 766.6*U.AA,
            stardust_code = 317
        ),
        Filter(
            instrument = 'HSC',
            telescope = 'Subaru',
            band = 'y',
            spec_wavelength = 9762*U.AA,
            spec_width = 787.5*U.AA,
            stardust_code = 318
        ),
        Filter(
            instrument = 'VIRCAM',
            telescope = 'VISTA',
            band = 'Y',
            spec_wavelength = 10190*U.AA,
            spec_width = 924.5*U.AA,
            stardust_code = 256
        ),
        Filter(
            instrument = 'VIRCAM',
            telescope = 'VISTA',
            band = 'J',
            spec_wavelength = 12460*U.AA,
            spec_width = 1718*U.AA,
            stardust_code = 257
        ),
        Filter(
            instrument = 'VIRCAM',
            telescope = 'VISTA',
            band = 'H',
            spec_wavelength = 16310*U.AA,
            spec_width = 2903*U.AA,
            stardust_code = 258
        ),
        Filter(
            instrument = 'VIRCAM',
            telescope = 'VISTA',
            band = 'Ks',
            spec_wavelength = 21400*U.AA,
            spec_width = 3076*U.AA,
            stardust_code = 259
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'IB427',
            spec_wavelength = 4255*U.AA,
            spec_width = 207.3*U.AA,
            stardust_code = 181
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'IB464',
            spec_wavelength = 4632*U.AA,
            spec_width = 218.2*U.AA,
            stardust_code = 183
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'IA484',
            spec_wavelength = 4845*U.AA,
            spec_width = 229.2*U.AA,
            stardust_code = 184
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'IB505',
            spec_wavelength = 5060*U.AA,
            spec_width = 231.6*U.AA,
            stardust_code = 185
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'IA527',
            spec_wavelength = 5258*U.AA,
            spec_width = 242.9*U.AA,
            stardust_code = 186
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'IB574',
            spec_wavelength = 5761*U.AA,
            spec_width = 272.9*U.AA,
            stardust_code = 188
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'IA624',
            spec_wavelength = 6228*U.AA,
            spec_width = 300.4*U.AA,
            stardust_code = 190
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'IA679',
            spec_wavelength = 6777*U.AA,
            spec_width = 336.3*U.AA,
            stardust_code = 192
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'IB709',
            spec_wavelength = 7069*U.AA,
            spec_width = 316.3*U.AA,
            stardust_code = 193
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'IA738',
            spec_wavelength = 7680*U.AA,
            spec_width = 364.8*U.AA,
            stardust_code = 194
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'IA767',
            spec_wavelength = 7357*U.AA,
            spec_width = 323.5*U.AA,
            stardust_code = 195
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'IB827',
            spec_wavelength = 8240*U.AA,
            spec_width = 343*U.AA,
            stardust_code = 197
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'NB711',
            spec_wavelength = 7120*U.AA,
            spec_width = 72.57*U.AA,
            stardust_code = 322
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'NB816',
            spec_wavelength = 8149*U.AA,
            spec_width = 119.8*U.AA,
            stardust_code = 319
        ),
        Filter(
            instrument = 'VIRCAM',
            telescope = 'VISTA',
            band = 'NB118',
            spec_wavelength = 11910*U.AA,
            spec_width = 112.2*U.AA,
            stardust_code = 321
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'B',
            spec_wavelength = 4420*U.AA,
            spec_width = 893.3*U.AA,
            stardust_code = 114
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'V',
            spec_wavelength = 5434*U.AA,
            spec_width = 954.6*U.AA,
            stardust_code = 115
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'r+',
            spec_wavelength = 6205*U.AA,
            spec_width = 1376*U.AA,
            stardust_code = 116
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'i+',
            spec_wavelength = 7604*U.AA,
            spec_width = 1497*U.AA,
            stardust_code = 117
        ),
        Filter(
            instrument = 'Suprime-Cam',
            telescope = 'Subaru',
            band = 'z++',
            spec_wavelength = 9070*U.AA,
            spec_width = 1335*U.AA,
            stardust_code = 118
        ),
        Filter(
            instrument = 'IRAC',
            telescope = 'Spitzer',
            band = 'ch1',
            spec_wavelength = 35130*U.AA,
            spec_width = 7445*U.AA,
            stardust_code = 18
        ),
        Filter(
            instrument = 'IRAC',
            telescope = 'Spitzer',
            band = 'ch2',
            spec_wavelength = 44430*U.AA,
            spec_width = 10120*U.AA,
            stardust_code = 19
        ),
        Filter(
            instrument = 'IRAC',
            telescope = 'Spitzer',
            band = 'ch3',
            spec_wavelength = 56480*U.AA,
            spec_width = 14070*U.AA,
            stardust_code = 20
        ),
        Filter(
            instrument = 'IRAC',
            telescope = 'Spitzer',
            band = 'ch4',
            spec_wavelength = 76180*U.AA,
            spec_width = 28790*U.AA,
            stardust_code = 21
        ),
        Filter(
            instrument = 'GALEX',
            telescope = 'GALEX',
            band = 'FUV',
            spec_wavelength = 1543*U.AA,
            spec_width = 227.8*U.AA,
            stardust_code = 120
        ),
        Filter(
            instrument = 'GALEX',
            telescope = 'GALEX',
            band = 'NUV',
            spec_wavelength = 2277*U.AA,
            spec_width = 790.9*U.AA,
            stardust_code = 121
        )
    ]
    return stardust_dotspec_map

if __name__ == '__main__':
    a = get_cosmos_dotspec_map()
    print(len(a))