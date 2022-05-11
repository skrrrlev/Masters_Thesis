
from enum import Enum

class PATH(Enum):
    ROOT = '/Volumes/Git/Masters_Thesis/'
    FIGURES = ROOT + 'Figures/'
    DATA = ROOT + 'Data/'
    DATAMAPS = DATA + 'Maps/'
    PHOTOMETRY = DATA + 'Photometry/'
    SCRIPTS = ROOT + 'Scripts/'
    OUTPUT = ROOT + 'Output/'
    OUTPUTMAPS = OUTPUT + 'Maps/'

if __name__ == '__main__':
    print(PATH.ROOT.value)
    