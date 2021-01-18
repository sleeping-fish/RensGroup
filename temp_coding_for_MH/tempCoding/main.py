# @Author:MH, LXY
# -*- coding = 'utf-8' -*-
# @Time:2020/11/20 08:39
# @File: main.py
# @Software: PyCharm


'''
In this programm, pdb files are read in to generate an instance of the protein data structure.
The direction of a residue can be further calculated.
'''

import glob
import re
from positioner import Protein


if __name__ == '__main__':
    pdbFiles = glob.glob(r"D:\pythonProject\Python_code\temp_coding_for_MH\tempCoding\PdbSet\*.pdb")
    XlsxDir = r'D:\pythonProject\Python_code\temp_coding_for_MH\tempCoding\XlsxDir\\'
    TxtDir = r'D:\pythonProject\Python_code\temp_coding_for_MH\tempCoding\TxtSet'
    SavePttern = re.compile(r'(.*?).pdb')

    for filename in pdbFiles:
        Name = re.findall(SavePttern, filename)
        SaveNames = Name[0][-3:]+'.txt'
        myStructure = Protein(filename, XlsxDir, TxtDir, "myProtein")
        myStructure.toStructure()
        myStructure.getResidue(SaveNames)

    print('Finished!')
