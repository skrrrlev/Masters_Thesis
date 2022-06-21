
from sys import argv
from os.path import isdir

def setup():
    root = argv[1]
    if not isdir(root):
        raise ValueError('Root must be a directory.')
    return root


def main():
    pass




if __name__ == '__main__':
    main()