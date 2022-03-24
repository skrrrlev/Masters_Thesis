# imports
import stardust
import os
import inspect
from re import findall

class CustomStardustFilters:

    rep = os.path.abspath(inspect.getfile(stardust.main.ctf)).strip('main.py') + 'filters/'
    filters = rep + 'filters.txt'
    info = rep + 'filters.info'


    @classmethod
    def get_number_of_filters(cls):
        '''Read number of filters in stardust from the filters.info file.'''
        with open(CustomStardustFilters.info, 'r') as fi:
            contents = fi.read()
        lines = contents.split('\n')
        return len(lines)

    @classmethod
    def is_filter_in_info_file(cls, filter: str) -> bool:
        '''Check if a filter is already in the filters.info file.'''
        with open(CustomStardustFilters.info, 'r') as fi:
            contents = fi.read()
        lines = contents.split('\n')
        flag = False
        for line in lines:
            if filter in line:
                flag = True
        return flag

    @classmethod
    def add_file(cls, file: str):
        '''Add filters from a file to the filters.info and filters.txt files.'''
        with open(file, 'r') as f:
            contents = f.read()
            lines = contents.split('\n')
        
        filters = []
        current_filter = []
        lines_in_filter = 0
        for line in lines:
            if not current_filter:
                lines_in_filter = int(findall(r'[\s+]?(\d+)[\s\S]*',line)[0])
            current_filter.append(line.strip())
            if len(current_filter)-1 == lines_in_filter:
                filters.append(current_filter)
                current_filter = []
        
        print(f'Found {len(filters)} filters in {file}')
        
        number = CustomStardustFilters.get_number_of_filters() + 1
        for filter in filters:
            info = filter[0]
            if not CustomStardustFilters.is_filter_in_info_file(info):
                print(f'Added filter no. {number} to stardust filter list: {info}')
                with open(CustomStardustFilters.info, 'a+') as fa:
                    fa.write(f'\n{number} {info}')
                with open(CustomStardustFilters.filters, 'a+') as fa:
                    fa.write(f'\n{info}')
                    for line in filter[1:]:
                        fa.write(f'\n{line}')
                number += 1
            else:
                print(f'Filter was already present in stardust filter list: {info}')




if __name__ == '__main__':
    c = CustomStardustFilters.get_number_of_filters()
    CustomStardustFilters.add_file('Data/COSMOS/extra_filters.txt')