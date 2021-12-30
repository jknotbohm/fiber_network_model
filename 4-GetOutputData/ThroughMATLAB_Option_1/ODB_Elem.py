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
    NOI=inst[(inst.keys()[i])].elementSets[NS].elements    
    #break
  except:
    pass

E_Label = NOI[0].label
E_Nodes = NOI[0].connectivity
E_inst  = NOI[0].instanceName
# Get Element Length
N1 = inst[E_inst].nodes[E_Nodes[0]-1]
N2 = inst[E_inst].nodes[E_Nodes[1]-1]
L=0
L_vec = N1.coordinates-N2.coordinates
for i in L_vec:
  L = L + i*i
L = sqrt(L)


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
FO_line=[];
for j in range(0,len(FO)):
  FO_full = FinalFrame.fieldOutputs[FO[j]]
  try:
    FO_line.append(FO_full.values[-1].data)
  except:
    FO_line.append(array('f',[FO_full.values[-1].data]))
FO_set.append(FO_line)


# Write Output To .txt file
txtName=odbName[0:-4]+'.txt'
print(txtName)
txtFile=open(txtName,'w')
txtFile.write("%10s,"%('Elem-Num') )
txtFile.write("%10s,"%('L') )
for i in range(0,L_FO):
  txtFile.write("%10s,"%(FO_labels[i]))
txtFile.write("\n")


txtFile.write("%10d,"%(E_Label))
txtFile.write("%10f,"%(L))
for j in range (0,len(FO_set[0])):
  for k in range(0,len(FO_set[0][j])):
    txtFile.write("%10f,"%(FO_set[0][j][k]))
txtFile.write("\n")

#Clean Up
print  #Just add a blank line
myOdb.close()
txtFile.close()