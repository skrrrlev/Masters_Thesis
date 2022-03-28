# imports
from typing import NoReturn
import sys
from os.path import isdir, isfile
from astropy import units as U
import numpy as np
import string
alpha_dict = dict(enumerate(string.ascii_lowercase))


# MTLib
from MTLib.data import xmlp


def usage() -> NoReturn:
    docs = """
    python table_from_xml.py <xml_file> [<output_path> [<output_name>]]
        <xml_file>    : path to xml file
            
        <output_path> : path to output
            If output_path is not given, the path will defeault to directory of the xml file.
        <output_name> : name of output
            If output_name is not given, the name will default to the name of the xml file.
    """
    print(docs)
    exit()

def setup() -> "tuple[str,str,str]":
    number_of_args = len(sys.argv) - 1
    if not number_of_args or number_of_args > 3:
        usage()
    print(number_of_args)

    xml_file: str = sys.argv[1].replace('\\','/')
    if not isfile(xml_file):
        usage()
    print(xml_file)

    if number_of_args > 1:
        output_path: str = sys.argv[2]
        if not isdir(output_path):
            print(f'{output_path} is not a directory')
            usage()
    else:
        output_path = '/'.join(xml_file.split('/')[:-1])
    print(output_path)

    if number_of_args == 3:
        output_name: str = sys.argv[3]
    else:
        output_name: str = xml_file.split('/')[-1].strip('.xml')
    print(output_name)
    
    return xml_file, output_path, output_name

to = 'uJy'
def convert(unit,f,Δf):
    if unit == to:
        return f,Δf
    elif 'Jy' in unit:
        f = ((f*U.Unit(unit)).to(to)).value
        Δf = ((Δf*U.Unit(unit)).to(to)).value
        return f,Δf
    elif unit =='ABmag':
        df_dm = - ( 3631 * np.log(10) * np.exp(-f * np.log(10) / 2.5) ) / 2.5
        Δf = (np.sqrt((df_dm * Δf)**2)*U.Jy).to(to).value
        f = ((f*U.ABmag).to(to)).value
        return f,Δf
    else:
        raise NotImplementedError(f'Invalid unit: "{unit}"')


def get_latex_unit(unit:str):
    if unit == 'mJy':
        return '\\si{\\milli\\jansky}'
    elif unit == 'uJy':
        return '\\si{\\micro\\jansky}'
    elif unit == 'ABmag':
        return '\\si{\\mag}'
    else:
        raise NotImplementedError(f'Invalid unit: "{unit}"')

def main():
    xml_file, output_path, output_name = setup()

    targets = xmlp.Parser.load_targets(xml_file)
    sources = xmlp.Parser.load_sources(xml_file)

    id_to_name = {str(int(target['id'])):target['name'] for target in targets}

    target_list = []
    obs_list = []
    columns = {}
    source_number: dict[str,int] = {}
    for source in sources:
        n = len(source_number)
        source_number[source] = alpha_dict[n]
        for observation in sources[source]:
            name_of_observation = observation['name']
            target_of_observation = id_to_name[observation['target']]

            if not target_of_observation in target_list:
                target_list.append(target_of_observation)

            if not name_of_observation in obs_list:
                obs_list.append(name_of_observation)

            f,Δf = convert(observation['unit'],float(observation['flux']),float(observation['error']))
            if f == Δf:
                fl = f'{f:.1e}'.split('e')
                string = f'< {fl[0]}E${"+" if (int(fl[1])>=0) else "-"}${abs(int(fl[1]))} $^{alpha_dict[n]}$'
            else:
                fl = f'{f:.1e}'.split('e')
                Δfl = f'{Δf:.1e}'.split('e')
                string = f' {fl[0]}Ε${"+" if (int(fl[1])>=0) else "-"}${abs(int(fl[1]))} ± {Δfl[0]}Ε${"+" if (int(Δfl[1])>=0) else "-"}${abs(int(Δfl[1]))} $^{alpha_dict[n]}$'
            if not target_of_observation in columns:
                columns[target_of_observation] = {name_of_observation:string}
            else:
                columns[target_of_observation][name_of_observation] = string
    
    number_of_obs = len(obs_list)
    number_of_targets = len(target_list)
    table = '\\begin{ThreePartTable}\n\\begin{TableNotes}\n\\footnotesize\n'
    for source in source_number:
        table += f'\\item [${source_number[source]}$] {source}\n'
    table += '\\end{TableNotes}\n\\begin{footnotesize}\n\\begin{longtable}{@{}l'+"c"*number_of_targets+'@{}}\n'
    table += '\\caption{Flux densities in ' + f'${get_latex_unit(to)}$. '
    for i,source in enumerate(source_number):
        table += f'${source_number[source]}$:{source}'
        if i != len(source_number)-1:
            table += ', '
        else:
            table += '.}\n'
    table += '\\label{tab:cat}\\\\\n\\toprule\n'
    table += 'Observation & '
    number_of_cols = len(columns)
    for i,col in enumerate(columns):
        table += f'{col}'
        if i+1 != number_of_cols:
            table += ' & '
        else:
            table += ' \\\\* \\midrule\n'
    table += '\\endfirsthead\n%\n\\endhead\n%\n\\bottomrule\n\\endfoot\n\\bottomrule\n\\insertTableNotes\n%\n\\endlastfoot\n%\n'
    for i,obs in enumerate(obs_list):
        table += obs.replace('um','')
        for column in columns:
            try:
                table += f' & {columns[column][obs]}'
            except KeyError:
                table += ' & $\\cdots$'
        table += ' \\\\'
        if i == len(obs_list)-1:
            table += '* \\bottomrule'
        table += '\n'
    table += '\\end{longtable}\n\\end{footnotesize}\n\\end{ThreePartTable}'
    print(table)

if __name__ == '__main__':
    main()