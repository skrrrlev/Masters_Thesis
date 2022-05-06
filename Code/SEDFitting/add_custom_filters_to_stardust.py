"""
author: Ditlev Frickmann
email: reimer.frickmann@gmail.com

Before stardust is envoked in fitting, this script needs to be run, to add the extra filters defined in the given file.

"""
from MTLib.data.stardust import CustomStardustFilters as csf

def main():
    csf.add_file('Data/stardust/extra_filters.txt')

if __name__ == '__main__':
    main()