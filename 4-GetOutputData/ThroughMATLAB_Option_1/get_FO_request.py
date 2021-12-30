# This function will report nodal position and reaction forces 
# for specified nodeset

from __future__ import division
from odbAccess import *
import odbAccess
from sys import argv,exit
from math import fabs
import sys 
from sys import path
import os
import platform

# Variable Setup
part = 'PART-1-1'
odbName=sys.argv[1]

# Get Current Working Dir
baseDir = os.getcwd();

# Open Odb - Last Step - Last Frame
myOdb = openOdb(baseDir+'/'+odbName)
AllSteps=myOdb.steps
FinalStep=myOdb.steps[str(AllSteps.keys()[-1])]
FinalFrame=FinalStep.frames[-1]
AllKeys = FinalFrame.fieldOutputs.keys()
for i in range(0,len(AllKeys)):
  print AllKeys[i]
  

