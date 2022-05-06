from sys import setrecursionlimit, getrecursionlimit
from os.path import isdir
from os import listdir


class recursion_limit:
    '''Set the recursion limit for the current scope'''
    def __init__(self, scope_limit:int = 10000):
        self.out_of_scope_limit = getrecursionlimit()
        self.scope_limit = scope_limit

    def __enter__(self):
        setrecursionlimit(self.scope_limit)

    def __exit__(self, exc_type, exc_val, exc_tb):
        setrecursionlimit(self.out_of_scope_limit)


def __walker(root: str, extension='*'):
    root = root.replace('\\','/')
    root += '/'*(not root.endswith('/'))
    files = []
    for item in listdir(root):
        path = root + item
        if isdir(path):
            files += __walker(
                root=root+item,
                extension=extension
            )
        elif extension == '*':
            files.append(path)
        elif item.endswith(extension):
            files.append(path)
    return files

def walker(root: str, extension='*'):
    '''Specify a root directory and get all files, within the directory and all directories under the root, with a given extension.
    
    If extension='*', then all files are returned.'''
    with recursion_limit(10000):
        return __walker(root,extension)

if __name__ == '__main__':
    files = walker('Data/Maps/', '.reg')
    print(files,'\n')
    files = walker('Data/Maps/', '.fits')
    print(files)