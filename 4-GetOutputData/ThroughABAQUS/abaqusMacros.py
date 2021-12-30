# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__


def DirectoryOutput_El():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import os
    import glob, os
    Main_Dir = "F:/Stephen/2019-BurkelHelp/2019-07-22 Wall Post/SampleFiles/"
    dir_walk = os.walk(Main_Dir)
    subdirs = next(dir_walk)[1]
    for s in subdirs:
        Working_Dir = ''.join(Main_Dir+s+"/")
        print(Working_Dir)
        os.chdir(Working_Dir)
        for odb_list in glob.glob("*.odb"):
            odb_name = str(odb_list)
            print(odb_name)
            rpt_name = odb_name[:(len(odb_name)-4)]
            rpt_name =''.join(rpt_name+"_Elemental.rpt")
            file_loc = ''.join(Working_Dir+odb_name)
            o1 = session.openOdb(
                name=file_loc)
            odb = session.odbs[file_loc]
            step_temp = odb.steps[str(odb.steps.keys()[-1])]
            frames = step_temp.frames
            session.viewports['Viewport: 1'].setValues(displayedObject=o1)

            session.writeFieldReport(rpt_name, append=OFF, 
                sortItem='Element Label', odb=odb, step=len(odb.steps)-1, frame=len(step_temp.frames)-1, 
                outputPosition=INTEGRATION_POINT, variable=(('S', INTEGRATION_POINT, ((
                COMPONENT, 'S11'), (COMPONENT, 'S22'), (COMPONENT, 'S33'), )), ))

            
            session.odbs[file_loc].close(
            )
            
def FolderOutput_El():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import os
    import glob, os
    Working_Dir = "F:/PaperReview_NetProperties/2019-2-12 Bending Full/"
    print(Working_Dir)
    os.chdir(Working_Dir)
    for odb_list in glob.glob("*.odb"):
        odb_name = str(odb_list)
        print(odb_name)
        rpt_name = odb_name[:(len(odb_name)-4)]
        rpt_name =''.join("Elemental_"+rpt_name+".rpt")
        file_loc = ''.join(Working_Dir+odb_name)
        o1 = session.openOdb(
            name=file_loc)
        odb = session.odbs[file_loc]
        step_temp = odb.steps[str(odb.steps.keys()[-1])]
        frames = step_temp.frames
        session.viewports['Viewport: 1'].setValues(displayedObject=o1)

        session.writeFieldReport(rpt_name, append=OFF, 
            sortItem='Element Label', odb=odb, step=len(odb.steps)-1, frame=len(step_temp.frames)-1, 
            outputPosition=INTEGRATION_POINT, variable=(('SE', INTEGRATION_POINT),('SK', INTEGRATION_POINT)),)
        
        session.odbs[file_loc].close(
        )

def DirectoryOutput_Node():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import os
    import glob, os
    
    Main_Dir = "F:/PaperReview_NetProperties/Density/"
    dir_walk = os.walk(Main_Dir)
    subdirs = next(dir_walk)[1]
    for s in subdirs:
        Working_Dir = ''.join(Main_Dir+s+"/Uniaxial Tension Tests/")
        print(Working_Dir)
        os.chdir(Working_Dir)
        for odb_list in glob.glob("*.odb"):
            odb_name = str(odb_list)
            print(odb_name)
            rpt_name = odb_name[:(len(odb_name)-4)]
            rpt_name =''.join("DISP_"+rpt_name+".rpt")
            file_loc = ''.join(Working_Dir+odb_name)
            o1 = session.openOdb(
                name=file_loc)
            odb = session.odbs[file_loc]
            step_temp = odb.steps[str(odb.steps.keys()[-1])]
            frames = step_temp.frames
            session.viewports['Viewport: 1'].setValues(displayedObject=o1)

            session.writeFieldReport(rpt_name, append=OFF, 
                sortItem='Node Label', odb=odb, step=len(odb.steps)-1, frame=len(step_temp.frames)-1, 
                outputPosition=NODAL, variable=(('U', NODAL),('UR', NODAL)),)
            
            session.odbs[file_loc].close(
                )
                
                

def FolderOutput_Node():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import os
    import glob, os
    
    Working_Dir = "F:/PaperReview_NetProperties/2019-2-12 Bending Full/"
    print(Working_Dir)
    os.chdir(Working_Dir)
    for odb_list in glob.glob("*.odb"):
        odb_name = str(odb_list)
        print(odb_name)
        rpt_name = odb_name[:(len(odb_name)-4)]
        rpt_name =''.join("DISP_"+rpt_name+".rpt")
        file_loc = ''.join(Working_Dir+odb_name)
        o1 = session.openOdb(
            name=file_loc)
        odb = session.odbs[file_loc]
        step_temp = odb.steps[str(odb.steps.keys()[-1])]
        frames = step_temp.frames
        session.viewports['Viewport: 1'].setValues(displayedObject=o1)

        session.writeFieldReport(rpt_name, append=OFF, 
            sortItem='Node Label', odb=odb, step=len(odb.steps)-1, frame=len(step_temp.frames)-1, 
            outputPosition=NODAL, variable=(('U', NODAL),('UR', NODAL)),)
        
        session.odbs[file_loc].close(
            )
            