# imports
from astropy.io import fits
from astropy import units as U
from typing import Union

# relative imports
from ..dataclasses import Observation, Filter

id_map: "dict[int,Union[int,None]]"= {
    926907: 911594,
    656648: 877661,
    82551: 832735,
    572380: 777178,
    659239: 884363,
    30861198: None,
    135815: 1452625,
    674217: 970128,
    673258: 970358,
    10131077: 972923,
}
'''Link the IDs from the deblended fits file to the IDs of the spec files.'''

def get_data(unit: U.Quantity=U.mJy) -> "list[Observation]":
    '''extract data from the super-deblended fits file'''

    observations: list[Observation] = []

    with fits.open('Data/COSMOS/Super_deblended_photometry_and_A3COSMOS.fits') as f:
        ir: fits.hdu.table.BinTableHDU = f[1]
        ir_dict = {line[0]:{'RA':line[1],'DEC':line[2]} for line in ir.data}
        ir_data = ir.data

    keys = ['SPITZER/MIPS/24','Herschel/PACS/100','Herschel/PACS/160','Herschel/SPIRE/250','JCMT/SCUBA-2/850','ALMA/870','ALMA/1000','ALMA/1250','ALMA/3000','VLA/3GHz','VLA/1.5GHz','Meerkat/1.3GHz']

    λ = {
        'SPITZER/MIPS/24': 24*U.um,
        'Herschel/PACS/100': 100*U.um,
        'Herschel/PACS/160': 160*U.um,
        'Herschel/SPIRE/250': 250*U.um,
        'JCMT/SCUBA-2/850': 850*U.um,
        'ALMA/870': 870*U.um,  
        'ALMA/1000': 1000*U.um,
        'ALMA/1250': 1250*U.um,
        'ALMA/3000': 3000*U.um,
        'VLA/3GHz': 10*U.cm,
        'VLA/1.5GHz': 20*U.cm,
        'Meerkat/1.3GHz': 23*U.cm,
    }
    code: dict[str,Filter] = {
        'SPITZER/MIPS/24': Filter('MIPS','SPITZER','24',850*U.um,0*U.m,325),
        'Herschel/PACS/100': Filter('PACS','Herschel','100',100*U.um,0*U.m,329), # Herschel/PACS/100
        'Herschel/PACS/160': Filter('PACS','Herschel','160',160*U.um,0*U.m,330), # Herschel/PACS/100
        'Herschel/SPIRE/250': Filter('SPIRE','Herschel','250',250*U.um,0*U.m,331), # Herschel/SPIRE/250um 331
        'JCMT/SCUBA-2/850': Filter('SCUBA-2','JCMT','850um',850*U.um,0*U.m,324),
        'ALMA/870': None,  
        'ALMA/1000': None,
        'ALMA/1250': None,
        'ALMA/3000': None,
        'VLA/3GHz': None,
        'VLA/1.5GHz': None,
        'Meerkat/1.3GHz': None,
    }

    for line in ir_data:
        
        idx = id_map[line[0]]
        
        if isinstance(idx,tuple):
            idx = idx[0]
        if idx is None:
            continue
        
        obs = {
            'SPITZER/MIPS/24': line[15] * U.uJy, # MIPS Spitzer
            'Herschel/PACS/100': line[17] * U.mJy, # Herschel PACS
            'Herschel/PACS/160': line[19] * U.mJy, # Herschel PACS
            'Herschel/SPIRE/250': line[21] * U.mJy, # spire 250 μm
            # SCUBA2
            'JCMT/SCUBA-2/850': line[23] * U.mJy,
            # ALMA
            'ALMA/870': line[25] * U.mJy,  
            'ALMA/1000': line[27] * U.mJy,
            'ALMA/1250': line[29] * U.mJy,
            'ALMA/3000': line[31] * U.mJy,
            # VLA
            'VLA/3GHz': line[33] * U.uJy,
            'VLA/1.5GHz': line[35] * U.uJy,
            # Meerkat
            'Meerkat/1.3GHz': line[37] * U.uJy,
        }

        err = {
            'SPITZER/MIPS/24': line[16] * U.uJy,
            'Herschel/PACS/100': line[18] * U.mJy,
            'Herschel/PACS/160': line[20] * U.mJy,
            'Herschel/SPIRE/250': line[22] * U.mJy,
            'JCMT/SCUBA-2/850': line[24] * U.mJy,
            'ALMA/870': line[26] * U.mJy,
            'ALMA/1000': line[28] * U.mJy,
            'ALMA/1250': line[30] * U.mJy,
            'ALMA/3000': line[32] * U.mJy,
            'VLA/3GHz': line[34] * U.uJy,
            'VLA/1.5GHz': line[36] * U.uJy,
            'Meerkat/1.3GHz': line[37] * U.uJy,
        }
        for key in keys:
            if obs[key].value == -99:
                continue
            if code[key] is None:
                new_obs = Observation(
                    target_id = idx,
                    flux = obs[key].to(unit).value,
                    flux_error= err[key].to(unit).value,
                    unit = unit,
                    wavelength = λ[key].to(U.um)
                )
                new_obs.set_name(key)
            else:
                new_obs = Observation(
                    target_id = idx,
                    flux = obs[key].to(unit).value,
                    flux_error= err[key].to(unit).value,
                    unit = unit,
                    filter = code[key]
                )
            observations.append(new_obs)
    
    return observations