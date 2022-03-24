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

    observations: list[Observation] = []

    with fits.open('Data/COSMOS/Super_deblended_photometry_and_A3COSMOS.fits') as f:
        ir: fits.hdu.table.BinTableHDU = f[1]
        ir_dict = {line[0]:{'RA':line[1],'DEC':line[2]} for line in ir.data}
        ir_data = ir.data

    keys = ['f24','f100','f160','f250','f850','f870','f1000','f1250','f3000','f10cm','f20cm','f23cm']

    λ = {
        'f24': 24*U.um,
        'f100': 100*U.um,
        'f160': 160*U.um,
        'f250': 250*U.um,
        'f850': 850*U.um,
        'f870': 870*U.um,  
        'f1000': 1000*U.um,
        'f1250': 1250*U.um,
        'f3000': 3000*U.um,
        'f10cm': 10*U.cm,
        'f20cm': 20*U.cm,
        'f23cm': 23*U.cm,
    }
    code = {
        'f24': None,
        'f100': None,
        'f160': None,
        'f250': None,
        'f850': Filter('SCUBA-2','JCMT','850um',850*U.um,0*U.m,324),
        'f870': None,  
        'f1000': None,
        'f1250': None,
        'f3000': None,
        'f10cm': None,
        'f20cm': None,
        'f23cm': None,
    }

    for line in ir_data:
        
        idx = id_map[line[0]]
        
        if isinstance(idx,tuple):
            idx = idx[0]
        if idx is None:
            continue
        
        obs = {
            'f24': line[15] * U.uJy, # MIPS Spitzer
            'f100': line[17] * U.mJy, # Herschel PACS
            'f160': line[19] * U.mJy, # Herschel PACS
            'f250': line[21] * U.mJy, # spire 250 μm
            # SCUBA2
            'f850': line[23] * U.mJy,
            # ALMA
            'f870': line[25] * U.mJy,  
            'f1000': line[27] * U.mJy,
            'f1250': line[29] * U.mJy,
            'f3000': line[31] * U.mJy,
            # VLA
            'f10cm': line[33] * U.uJy,
            'f20cm': line[35] * U.uJy,
            # Meerkat
            'f23cm': line[37] * U.uJy,
        }

        err = {
            'f24': line[16] * U.uJy,
            'f100': line[18] * U.mJy,
            'f160': line[20] * U.mJy,
            'f250': line[22] * U.mJy,
            'f850': line[24] * U.mJy,
            'f870': line[26] * U.mJy,
            'f1000': line[28] * U.mJy,
            'f1250': line[30] * U.mJy,
            'f3000': line[32] * U.mJy,
            'f10cm': line[34] * U.uJy,
            'f20cm': line[36] * U.uJy,
            'f23cm': line[37] * U.uJy,
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