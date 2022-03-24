
import xml.etree.ElementTree as ET
from typing import Union

class Parser:
    '''
    Parser that can be used to parse the data structures created for data in this thesis.

    The `load_targets` and `load_observations` methods are static methods that take an xml file as input.

    The Parser can be initialised to automatically add the targets and observations to the instance.  
    
    '''

    def __init__(self, file: str):
        self.targets = self.load_targets(file)
        self.observations = self.load_observations(file)

    @staticmethod
    def load_targets(file: str):
        tree = ET.parse(file)
        root = tree.getroot()
        targets_node = root[0]
        targets: "list[dict[str,Union[str,float]]]" = []

        for target in targets_node:
            temp = {}
            for item in target:
                try:
                    temp[item.tag] = float(item.text)
                except ValueError:
                    temp[item.tag] = item.text
            targets.append(temp)

        return targets

    @staticmethod
    def load_observations(file: str):
        tree = ET.parse(file)
        root = tree.getroot()
        sources = root[1]
        observations: "list[dict[str,Union[str,float]]]" = []

        for source in sources:
            for observation in source:
                temp = {}
                for key in observation.attrib:
                    temp[key] = observation.attrib[key]
                for item in observation:
                    try:
                        temp[item.tag] = float(item.text)
                    except ValueError:
                        temp[item.tag] = item.text
                observations.append(temp)

        return observations

    def get_targets(self):
        return self.targets

    def get_observations(self):
        return self.observations


if __name__ == '__main__':
    test = Parser('Data/XMM/XMM.xml')
    print('targets')
    print(test.get_observations())
    print('observations')
    print(Parser.load_observations('Data/XMM/XMM.xml'))