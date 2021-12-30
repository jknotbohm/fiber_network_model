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


#Creating a for loop to iterate through all frames in the step
for x in odb.steps[step1.name].frames:


####Setting temporary array to empty
	temp = []


#####Reading stress and strain data from the model 
	odbSelectResults = x.fieldOutputs['U']
	# odbSelectResults2 = x.fieldOutputs['U']

	field1 = odbSelectResults
	# field2 = odbSelectResults2

####Storing Stress and strain values for the current frame

	# stress values
	for s in field1.values:
		info = []
		temp.append(s.data[1]) # the values in s.data are organized as follows: [S11,S22,S33,S12,S13,S23]
	
	# strain values
	# for e in field2.values:
	#	temp.append(e.data[0]) # the values here are organized similarly to the stress values. This grabs E33

	DAT.append(temp)
	
    		
####Writing to a .csv file
with open(results_name+'.txt', 'w') as f:
	np.savetxt(f,DAT,delimiter='	')	


#Close the odb
odb.close()