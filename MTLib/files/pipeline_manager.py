
import configparser
from os import makedirs, remove
from os.path import isdir, isfile
from shutil import rmtree, copyfile

from . import walker

class NotInitializedError(Exception):
    '''Raised when you try to get the output of a .ini that has not been intialized.'''

class MapPipelineManager:

    def __init__(self, ini_file:str) -> None:
        self.ini = ini_file
        self.config = configparser.ConfigParser()
        self.config.read(self.ini)

    @staticmethod
    def read_input(argv) -> "list[str]":
        root: str = argv[1]
        if isdir(root):
            return walker(root,extension='.ini')
        elif isfile(root) and root.endswith('.ini'):
            return [root]
        else:
            raise ValueError('Root must be a directory or an .ini file.')

    def __enter__(self):
        return self
    
    def __exit__(self ,type, value, traceback):
        with open(self.ini,'w') as configfile:
            self.config.write(configfile)

    @staticmethod
    def initialize_from_file(ini_file) -> "list[str]":
        in_config = configparser.ConfigParser()
        in_config.read(ini_file)
        number_of_input = int(in_config['INPUT']['items'])

        new_files: "list[str]" = []

        for i in range(number_of_input):
            out_config = configparser.ConfigParser()
            
            # Load tab and get parameters
            try:
                tab = in_config['INPUT'][f'item{i}']
            except Exception as e:
                raise e(f'Expected n={number_of_input} number of input items in {ini_file}, but could not find item{i}.')
            try:
                input_source_file = in_config[tab]['source_file']
                input_region_file = in_config[tab]['region_file']
                output_path = in_config[tab]['output_path']
            except Exception as e:
                raise e(f'Could not initialize {tab}. It must contain the keywords "source_file", "region_file" and "output_path".')
            
            # Initialize output directory
            if isdir(output_path):
                rmtree(output_path)
            makedirs(output_path)
            output_ini = f'{output_path}/{tab}.ini'

            # Transfer the default values
            out_config['DEFAULT'] = in_config['DEFAULT']
            out_config['DEFAULT']['name'] = tab
            out_config['DEFAULT']['path'] = output_path
            
            # Move fits file and region file
            out_config['DEFAULT']['source_file'] = f'{output_path}/{tab}.fits'
            copyfile(input_source_file,out_config['DEFAULT']['source_file'])
            out_config['DEFAULT']['region_file'] = f'{output_path}/{tab}.reg'
            copyfile(input_region_file,out_config['DEFAULT']['region_file'])

            efm = 'exclude_from_mask'
            if efm in in_config[tab]:
                out_config['DEFAULT'][efm] = in_config[tab][efm]
            else:
                out_config['DEFAULT'][efm] = ''

            # Create the output ini file from the ConfigParser
            with open(output_ini, 'w') as f:
                out_config.write(f)
            new_files.append(output_ini)
        return new_files

    @staticmethod
    def get_output_files(ini_files: "list[str]"):
        '''For a list of input .ini files, get the corresponding output .ini files.'''
        files: "list[str]" = []
        for ini_file in ini_files:
            config = configparser.ConfigParser()
            config.read(ini_file)
            number_of_input = int(config['INPUT']['items'])
            for i in range(number_of_input):
                try:
                    tab = config['INPUT'][f'item{i}']
                except Exception as e:
                    raise e(f'Expected n={number_of_input} number of input items in {ini_file}, but could not find item{i}.')
                file = f'{config[tab]["output_path"]}/{tab}.ini'
                if not isfile(file):
                    raise NotInitializedError(f'Could not get output of {tab} in "{ini_file}", because the .ini has not been initialized with the "initialize_from_file" method. ("{file}" does not exist)')
                files.append(file)

        return files



    def get_output_tabs(self) -> "list[str]":
        '''Get the keys for the output tabs in the .ini file'''
        number_of_output = int(self.config['INPUT']['items'])
        return [self.config['OUTPUT'][f'item{i}'] for i in range(number_of_output)]

    def add_item_to(self, tab:str, item:str, value:str) -> None:
        '''Add <item> with <value> to <tab>.'''
        self.config[tab][item] = value

    def get_item_from(self, tab:str, item:str) -> str:
        '''Get <item> from <tab>.'''
        return self.config[tab][item]

    def contains(self, tab:str, item:str) -> bool:
        '''Returns whether <tab> contains <item>.'''
        return item in self.config[tab]

    def add_tab(self, tab):
        '''Add a new <tab>.'''
        self.config[tab] = {}

    def get_tab(self, tab:str):
        return self.config[tab]

if __name__ == '__main__':
    with MapPipelineManager('Data/Maps/HST/source14/source14.ini') as man:
        man.initialize()
        for key in man.config:
            print(key)

        output = man.get_output_tabs()
        print(output)
        man.add_item_to(output[0],'ayy','lmao')
        
