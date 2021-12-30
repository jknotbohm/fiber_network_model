import sys
from odbAccess import *
from abaqusConstants import*
from types import IntType
import numpy as np

# *.odb means every .odb file in the directory. Add directory path to the beginning: 'C:\path\*.odb'
odb_name='lattice_fiber_2020_07_05_trial3'

# desired name for the results file
results_name=odb_name+'test1'

#Initializing Arrays
DAT = []
info = []


#Opening the odb
odb = openOdb(odb_name+'.odb', readOnly=True)
assembly = odb.rootAssembly
instance = assembly.instances.keys()[0]

#Extracting Step 1, this analysis only had one step
step1 = odb.steps.values()[0]


for frame in odb.steps[step1.name].frames:
    print('%d,%f\n'%(frame.frameId,frame.frameValue))
    

#Close the odb
odb.close()
