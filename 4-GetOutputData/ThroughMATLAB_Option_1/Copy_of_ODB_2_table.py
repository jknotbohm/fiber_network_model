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
NS = sys.argv[2]
FOS = []
for i in range(3,len(sys.argv)):
 FOS.append(sys.argv[i])
print 'Model = ' , odbName
print 'Node Set = ', NS

# Get Current Working Dir
baseDir = os.getcwd();

# Open Odb - Last Step - Last Frame
myOdb = openOdb(baseDir+'/'+odbName)
AllSteps=myOdb.steps
FinalStep=myOdb.steps[str(AllSteps.keys()[-1])]
FinalFrame=FinalStep.frames[-1]
FO=[]
FOL=[]
for i in range(0,len(FOS)):
  FO.append(int(FOS[i])-1)
  FOL.append(FinalFrame.fieldOutputs.keys()[FO[i]])
print 'Field Outputs = ', FOL
inst=myOdb.rootAssembly.instances 


# Start Text File
txtName=odbName[0:-4]+'.txt'
print(txtName)
txtFile=open(txtName,'w')


# Get list of requested nodes
NS_list = []
NodeSet = inst[part].nodeSets[NS].nodes
for i in range(0,len(NodeSet)):
  NS_list.append(NodeSet[i].label)

# Get Get Nodal Positions and Displacements
NOI = inst[part].nodes
Node_Num=[]
Node_Coord=[]
for i in range(0,len(NOI)):
  if NOI[i].label in NS_list:
    Node_Num.append(NOI[i].label)
    Node_Coord.append(NOI[i].coordinates)
# Add To .txt File
txtFile.write("%10s,"%('NodeNum') )
txtFile.write("%10s,%10s,%10s,"%('X','Y','Z') )
txtFile.write("\n")
for i in range(0,len(Node_Num)):
  txtFile.write("%10s,"%(Node_Num[i]))
  for j in range(0,len(Node_Coord[i])):
    txtFile.write("%10s,"%(Node_Coord[i][j]))
  txtFile.write("\n")
txtFile.write("\n\n")  

# Repeat for ALL requested Feild Outputs
for f in range(0,len(FOL)):
  FOI=FinalFrame.fieldOutputs[FOL[f]]
  if len(FOI.componentLabels) == 0:
    FOI_Labels = [FOI.name];
  else:
    FOI_Labels = FOI.componentLabels
  txtFile.write("%10s,"%('NodeNum'))
  for i in range(0,len(FOI_Labels)):
    txtFile.write("%10s,"%(FOI_Labels[i]))
  txtFile.write("\n")
  for i in range(0,len(FOI.values)):
   if FOI.values[i].nodeLabel in NS_list:
    txtFile.write("%10s,"%(FOI.values[i].nodeLabel))
    try:
      len(FOI.values[i].data)
      FOI_Data = FOI.values[i].data
    except:
      FOI_Data = [FOI.values[i].data]
    for j in range(0,len(FOI_Data)):
      txtFile.write("%10s,"%(FOI_Data[j]))
    txtFile.write("\n")
  txtFile.write("\n\n")  


#Clean Up
print
myOdb.close()
txtFile.close()
