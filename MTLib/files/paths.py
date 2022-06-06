

def remove_path(file_path: str):
    '''Specify a path to a file and recieve the file without the path.'''
    file_path = file_path.replace('\\','/')
    return file_path[file_path.rfind('/')+1:]

def extract_path(file_path: str):
    '''Specify a path to a file and recieve the path to the directory of the file.'''
    file_path = file_path.replace('\\','/')
    return file_path[:file_path.rfind('/')+1]

def remove_extension(file_path: str):
    '''Specify a file or a path to a file and remove the extension of the file.'''
    return file_path[:file_path.rfind('.')]

def extract_extension(file_path: str):
    '''Specify a file or a path to a file and recieve the extension of the file.'''
    return file_path[file_path.rfind('.')+1:]

def extract_filename(file_path: str):
    '''Specify a path to a file and recieve the file name (path and extension is removed).'''
    return remove_extension(remove_path(file_path))

def extract_folder_name(folder: str):
    '''Specify a path to a folder and recieve the file name (path and extension is removed).'''
    folder = folder.replace('\\','/')
    if folder[-1] == '/':
        folder = folder[:-1]
    return folder[folder.rfind('/')+1:]

if __name__ == '__main__':
    file_path = 'Data/blabla/foo/bar/testfile.ext'
    print(f'file_path = {file_path}')
    print(f'Remove path:       {remove_path(file_path)}')
    print(f'Extract path:      {extract_path(file_path)}')
    print(f'Remove extension:  {remove_extension(file_path)}')
    print(f'Extract extension: {extract_extension(file_path)}')
    print(f'Extract file name: {extract_filename(file_path)}')
    print(f'Extract folder name: {extract_folder_name(extract_path(file_path))}')