#imports
from astropy import units as U
from os import listdir
from re import findall
import xml.etree.ElementTree as ET
from itertools import islice
from itertools import cycle

import sys

from C4S import Cataloguer
from C4S.dataclasses import ColumnType
from MTLib.data.dataclasses.observation import Observation

# private imports
from MTLib.data.stardust import CustomStardustFilters as csf
from MTLib.data.cosmos import SpecFile
from MTLib.data.cosmos import super_deblended as sd

cosmos_dir = 'Data/COSMOS/'
spec_input_dir = cosmos_dir+'COSMOS_uv_nir/'

split = 1
''' How many xml-files do you want? '''

cosmos_id_to_project_id = {
    877661:3,
    1024221:10,
    911594:29,
    832735:13,
    1452625:2,
    970128:20,
    884363:30,
    777178:24,
    970358:22,
    972923:19
}

def grouper(iterable,ngroups):
    groups = [[] for _ in range(ngroups)]
    for element, group in zip(iterable, cycle(groups)):
        group.append(element)
    return groups

def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

def main():

    """Load .spec files"""
    # Load .spec files and add them to fits file structure
    specfiles = [SpecFile(spec_input_dir + file, unit=U.mJy) for file in listdir(spec_input_dir) if file.endswith('.spec')]

    """Load RA and DEC and create target dictionary"""
    # open .txt file containing redshifts and RADEC, and read the contents into a string
    with open(cosmos_dir+'coord_cosmos_id2020.txt','r') as f:
        content = f.read()

    # Add them to lists.
    targets = {}
    for specfile in specfiles:
        targets[specfile.id] = {}

        ra_pattern = r'(\d+.\d+)\s+\d+.\d+\s+'+f'{specfile.id}'
        targets[specfile.id]['ra'] = float(findall(ra_pattern,content)[0])

        dec_pattern = r'\d+.\d+\s+(\d+.\d+)\s+'+f'{specfile.id}'
        targets[specfile.id]['dec'] = float(findall(dec_pattern,content)[0])

        targets[specfile.id]['z'] = specfile.z
        # targets[specfile.id]['z'] = specfile.get_zphot()

    target_ids_grouped = grouper(list(targets.keys()),split)

    for i,group in enumerate(target_ids_grouped):
        
        """Write beginning of xml file"""
        root = ET.Element("catalogue", name=f"CAT-COSMOS-{i+1}", label=f'COSMOS:photometry:{i+1}')
        targets_element = ET.SubElement(root,"targets")
        
        for target in targets:
            if target not in group:
                continue
            target_root = ET.SubElement(targets_element,"target")
            ET.SubElement(target_root,"name").text = f'{target}'
            ET.SubElement(target_root,"id").text = f'{cosmos_id_to_project_id[target]}'
            ET.SubElement(target_root,"RA").text = f'{targets[target]["ra"]}'
            ET.SubElement(target_root,"DEC").text = f'{targets[target]["dec"]}'
            ET.SubElement(target_root,"z").text = f'{targets[target]["z"]}'
        
        """Add observations from .spec files to catalogue"""
        sources_element = ET.SubElement(root,"data-sources")
        spec_source = ET.SubElement(sources_element,"source",name="Weaver:COSMOS2020")
        
        for specfile in specfiles:
            if specfile.id not in group:
                continue
            for observation in specfile.observations:
                # <observation target="24131301" name="ALMA-NOEMA/290um" unit="mJy" type="G">
                attrib = {"target":f"{cosmos_id_to_project_id[specfile.id]}","name":f"{observation.column.replace('_','/')}","unit":"mJy","type":"F"}
                obs_element = ET.SubElement(spec_source,"observation",attrib=attrib)
                ET.SubElement(obs_element,"flux").text = str(observation.flux)
                ET.SubElement(obs_element,"error").text = str(observation.flux_error)
                ET.SubElement(obs_element,"typeval").text = str(observation.filter.stardust_code)
        
        """ Load deblended data from fits file and add to xml"""
        deblended_source = ET.SubElement(sources_element,"source",name="Jin2018:super-deblended-photometry")
        a3cosmos_source = ET.SubElement(sources_element,"source",name="Liu:A3COSMOS-photometry")
        ir_data = sd.get_data(unit=U.mJy)
        for observation in ir_data:
            if observation.target_id not in group:
                continue
            if 'alma' in observation.column.lower() or 'herschel' in observation.column.lower():
                obs_source = a3cosmos_source
            else:
                obs_source = deblended_source
            if observation.type == ColumnType.FILTER:
                # <observation target="24131301" name="ALMA-NOEMA/290um" unit="mJy" type="F">
                attrib = {"target":f"{cosmos_id_to_project_id[observation.target_id]}","name":f"{observation.column.replace('_','/')}","unit":"mJy","type":"F"}
                obs_element = ET.SubElement(obs_source,"observation",attrib=attrib)
                ET.SubElement(obs_element,"flux").text = str(observation.flux)
                ET.SubElement(obs_element,"error").text = str(observation.flux_error)
                ET.SubElement(obs_element,"typeval").text = str(observation.filter.stardust_code)

            elif observation.type == ColumnType.EXTRA:
                # <observation target="24131301" name="ALMA-NOEMA/290um" unit="mJy" type="E">
                attrib = {"target":f"{cosmos_id_to_project_id[observation.target_id]}","name":f"{observation.column.replace('_','/')}","unit":"mJy","type":"E"}
                obs_element = ET.SubElement(obs_source,"observation",attrib=attrib)
                ET.SubElement(obs_element,"flux").text = str(observation.flux)
                ET.SubElement(obs_element,"error").text = str(observation.flux_error)
                ET.SubElement(obs_element,"typeval").text = str(observation.wavelength)
            else:
                raise ValueError('The type of the observation was not well defined')
        


        """Save file"""
        tree = ET.ElementTree(root)
        indent(root)
        tree.write(f'cosmos-{i+1}.xml',encoding="utf-8",xml_declaration=True)


if __name__ == '__main__':
    main()