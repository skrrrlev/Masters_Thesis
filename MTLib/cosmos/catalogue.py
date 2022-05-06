from astropy.io import fits
from os.path import dirname, abspath, isfile
import configparser
from typing import Union

class Catalogue:
    """
    Read data from the COSMOS catalogue used in this project.

    """
    def __init__(self) -> None:
        self.config = configparser.ConfigParser()
        self.inipath = dirname(abspath(__file__))
        self.config.read(self.inipath + '/properties.ini')

        with fits.open('Data/COSMOS/catalogue/COSMOS2020_CLASSIC_R1_v2.0.fits') as hdul:
            hdu: fits.hdu.table.BinTableHDU = hdul[1]
            self.data = hdu.data
            self.columns = hdu.columns

    def save_config_file(self):
        with open(self.inipath + '/properties.ini', 'w') as configfile:
            self.config.write(configfile)

    def add_to_ini(self, idxs: "Union[list[int], int]"):
        if isinstance(idxs,int):
            idxs = [idxs]

        objects_in_ini = list(self.config.keys())
        objects_not_in_ini = [idx for idx in idxs if str(idx) not in objects_in_ini]

        for idx in objects_not_in_ini:
            columns = self.data[idx-1]
            self.config[f'{columns[0]}'] = {}
            self.config[f'{columns[0]}']['lp_mass_med'] = str(columns[624])
            self.config[f'{columns[0]}']['lp_mass_best'] = str(columns[627])
            self.config[f'{columns[0]}']['ez_mass'] = str(columns[680])

        '''Regular
        624 name = 'lp_mass_med'; format = 'D'          log Stellar mass from BC03 best-fit template. median of the PDF 
        625 name = 'lp_mass_med_min68'; format = 'D'    lower limit, 68% confidence level
        626 name = 'lp_mass_med_max68'; format = 'D'    upper limit, 68% confidence level
        627 name = 'lp_mass_best'; format = 'D'         log Stellar mass from BC03 best-fit template. Taken at the minimum chi2
        680 name = 'ez_mass_p500'; format = 'D'; unit = 'solMass'   log(mass in Msun)
        '''        
            
    def get_stellar_masses(self, idxs: "Union[list[int],int]") -> "list[str]":
        if isinstance(idxs,int):
            idxs = [idxs]
        
        self.add_to_ini(idxs)
        return [self.config[str(idx)]['lp_mass_best'] for idx in idxs]

    def get_property(self,idxs: "Union[list[int],int]", property: str) -> "list[str]":
        if isinstance(idxs,int):
            idxs = [idxs]
        
        self.add_to_ini(idxs)
        return [float(self.config[str(idx)][property]) for idx in idxs]

    def print_columns(self):
        for i,col in enumerate(self.columns):
            print(i, col)


if __name__ == '__main__':
    cc = Catalogue()
    print(cc.get_stellar_masses([777178,832735,877661,884363,911594,970128,970358,972923,1024221,1452625]))
    cc.save_config_file()