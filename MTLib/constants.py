
from enum import Enum

class PATH(Enum):
    ROOT = __file__[:(__file__[:__file__.rfind('/')]).rfind('/')+1]
    FIGURES = ROOT + 'Figures/'
    DATA = ROOT + 'Data/'
    MAPS = DATA + 'Maps/'
    SCRIPTS = ROOT + 'Scripts/'

if __name__ == '__main__':
    print(PATH.ROOT.value)