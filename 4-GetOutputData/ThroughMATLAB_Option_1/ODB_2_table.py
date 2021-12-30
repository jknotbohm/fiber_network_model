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
from array  import array

# Variable Setup
odbName=sys.argv[1]
NS = sys.argv[2]
FO = []
for i in range(3,len(sys.argv)):
  FO.append(sys.argv[i])

print 'Model = ' , odbName
print 'NodeSet = ' , NS
print 'Field Outputs = ' , FO

# Get Current Working Dir
baseDir = os.getcwd();

# Open Odb - Last Step - Last Frame
myOdb = openOdb(baseDir+'/'+odbName)
AllSteps=myOdb.steps
FinalStep=myOdb.steps[str(AllSteps.keys()[-1])]
FinalFrame=FinalStep.frames[-1]


# Get Node Sets (Node Label & Coordinates)
inst=myOdb.rootAssembly.instances 
for i in range(0,len(inst.keys())):
  try:
    NOI=inst[(inst.keys()[i])].nodeSets[NS].nodes 
    #break
  except:
    pass

# Define Variables
Node_Num=[]; Node_Coord=[];
for i in range(0,len(NOI)):
  Node_Num.append(NOI[i].label)
  Node_Coord.append(NOI[i].coordinates)

# Get Desired Field Output
FO_full=[]; FO_temp=[]; FO_labels=[]
for i in range(0,len(FO)):
  FO_full = FinalFrame.fieldOutputs[FO[i]]
  FO_temp.append(FO_full.componentLabels)
for i in range(0,len(FO_temp)):
  for j in range(0,len(FO_temp[i])):
    FO_labels.append(FO_temp[i][j])
L_FO = len(FO_labels)
FO_set=[]; 
for i in Node_Num: 
  FO_line=[];
  for j in range(0,len(FO)):
    FO_full = FinalFrame.fieldOutputs[FO[j]]
    try:
      FO_line.append(FO_full.values[i-1].data)
    except:
      FO_line.append(array('f',[FO_full.values[i-1].data]))
  FO_set.append(FO_line)
# Check All Lengths Are The Same
L_labels=len(Node_Num)
L_Coord =len(Node_Coord)
L_Output=len(FO_set)

# Write Output To .txt file
txtName=odbName[0:-4]+'.txt'
print(txtName)
txtFile=open(txtName,'w')
txtFile.write("%10s,"%('NodeNum') )
txtFile.write("%10s,%10s,%10s,"%('X','Y','Z') )
for i in range(0,L_FO):
  txtFile.write("%10s,"%(FO_labels[i]))
txtFile.write("\n")
for i in range(0,L_labels):
  txtFile.write("%10d,"%(Node_Num[i]))
  txtFile.write("%10f,%10f,%10f,"%(Node_Coord[i][0],Node_Coord[i][1],Node_Coord[i][2]))
  for j in range (0,len(FO_set[i])):
    for k in range(0,len(FO_set[i][j])):
      txtFile.write("%.2E,"%(FO_set[i][j][k]))
  txtFile.write("\n")

#Clean Up
print  #Just add a blank line
myOdb.close()
txtFile.close()