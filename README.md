# READ ME

*Repository for fiber network simulations and modeling.*

*Notbohm Research Group, https://notbohm.ep.wisc.edu/*

This readme explains the Notbohm lab code to generate fiber networks and run 
simulations in Abaqus. This readme was originally written by Stephen Tyznik
with subsequent updates by other lab members. The initial version of the code 
was written by Peter Grimmer, with subsequent edits by Stephen Tyznik, Jacob 
Notbohm, and other lab members.

Names of Notbohm group members who contributed different sections of code
are written in the readme files. In general, the contributions are documented 
in the reserach manuscripts and theses written by these group members.

I highly recommend using the Abaqus and Matlab documentation whenever trying to
create something new or if you are confused on what anything does. Abaqus 6.14
documentation has a great keyword reference guide for all
commands we use. Matlab's online documentation is great as well.

# SECTION 1: NETWORK GENERATION

Original version written primarily by Peter Grimmer under guidance of Jacob Notbohm.

Optimized later by Stephen Tyznik. 

Subsequent updates and new capabilities added by Mainak Sarkar.

**The different methods are described along with documentation in a separate 
readme in the directory 1-NetworkGeneration.**

# SECTION 2: INP GENERATION

Each directory within 2-INPGeneration gives examples of inp generators
that apply different boundary conditions and employ different functionalities
in Abaqus. Each begins by loading a fiber network created during step 1.

**NOTE:** These are all just examples of what can be done with the fiber networks. It is extremely
important to use the ABAQUS documentation when trying to implement new features.

#### Tension+Bending
Written by Stephen Tyznik. Good intro to finite elements + Static Solver.

This is a simple displacement controlled boundary. 'Tension' pulls
the network along one axis. Allowing free movement and rotations along other
directions. ‘Bending’ is slightly more complicated and putes the network in a pure
bending case. Both use the STATIC solver.

#### Thermal Inclusion
Written by Stephen Tyznik. Example using continuum elements + Dynamic Solver.

A thin 3D network is flattened to a 2D network. A circular hole is
then cut in the center of the network. Continuum elements are designed to fill the hole
and are given a coefficient of thermal expansion. Thermal change is then applied and
shrinks the inclusion

#### Contact
Written by Stephen Tyznik. Example of contact problems.

Again a 2D network is flattened from a 3D network. The bottom boundary of
the fiber network is fixed and the upper is designated as a contact surface. Additionally
a rigid body is defined as a contact surface and pushed into the fiber network.

#### multi_step_INP_generators
Written by Mainak Sarkar. 

Newton-Raphson based quasistatic solvers used in Abaqus could not converge to force equilibrium in large deformation problems. Here we figured out an idea of 
imposing large deformation in multiple steps. This will reduce the convergence error at each incremental load step bringing down the overall error to a significant 
extent, if not completely eliminating it. A consideration of 20 steps was found to be sufficient for most purposes.

#### inp_gen_for lattice_network
Written by Mainak Sarkar.

`lattice_inp_main_script.m` and `inp_gen_lattice_network.m` are representative single-load-step INP file generation scripts written for uniaxial tensile boundary condition. 
This can guide users in scriptinng INP files related to lattice networks for other boundary conditions using multiple load steps. There is a supporting function file that finds the average fiber length 
in the network: `distance_finder.m`.

# SECTION 3: SUBMITTING JOBS

Two Examples for Submitting Jobs to Abaqus.

## GUI Based 
UI Interface to submit jobs through ABAQUS using a MATLAB window. Written by Stephen Tyznik.
`submit_abaqus_jobs`

##### Included Files
`Submit_Jobs_AllSubDir.m`

`Submit_Jobs_Folder.m`

##### Description
Both codes provided are essentially the same. With each you select all the folders with
INP files you want to run. Once done you are then able to use strings separated by
spaces to specify the specific files you run. An example of how to use this is to have
multiple windows open so you can run jobs in parallel.

The difference between the two codes is that `Submit_Jobs_Folder.m` will look only at
files in the immediate folder. Whereas `Submit_Jobs_AllSubDir.m` will look at the files in
the given directory AND all subdirectories.

##### Notes
**THIS IS NOT PERFECT.**
I wrote this as a way to improve efficiency but it has a lot of rough edges. The main two
being: 1) ‘completed’ does not take into account if there are errors, just that ABAQUS is
finished with the job. 2) Parallel computing cannot be done in the same window.
These were easy enough to work around that I never had time to fully fix them before I
left. Overall this code has still saved me a bunch of time so I wanted to include it.


## M-file Based 
A DOS based MATLAB-Abaqus interface to run an INP file in Abaqus and save the generated ODB file in the 
MATLAB's current directory.  
`run_inp_in_abaqus`

#### Included Files
`run_inp_in_abaqus_main_script.m`

`run_inp_abaqus.m`

##### Description
This technique is highly efficient and the scripts are concise, simple to implement in a loop of simulations.
It just runs Abaqus job from MATLAB through DOS, without even using the Abaqus User Interface.

# SECTION 4: GET OUTPUT DATA
A few options for getting ABAQUS .odb files into a useable format in matlab.

#### Through ABAQUS 
Written by Stephen Tyznik. Set of Abaqus macros to be run through the program.

Through ABAQUS is is an easy way to create rpt files directly from ABAQUS. Open up
the macro manager in the software and select either nodal or elemental outputs. In the
python code you can designate what outputs specifically you want. Like with submit
jobs, folder output will only run through a specific folders outputs while directory output
will go through all sub folders. You need to have this macro file saved in either the
Work or Home directory in Matlab.

**Disclaimer:** (Applies to "Through MATLAB" as well)
Again, this is a rough code I through together to help with my efficiency. It is by no
means perfect. I recommend using these as a template but learning what is going on
and either creating your own python code or abaqus macro to fit you specific need.
Rumor has it the Allen lab has code that does this as well. Might be worth the time to
ask them about it.

#### Through MATLAB (Option 1)
Written by Stephen Tyznik. Directly open the file through MATLAB and python.

Through MATLAB is a secondary option with two major benefits. First, you do not need
to use an Abaqus license to run this code. Second, you are able to specify a specific
node set rather than outputting for all nodes.

#### Through MATLAB (Option 2)
Written by Mainak Sarkar. Directly open the file through MATLAB and python.

This technique has all the benefits of previous method (Through MATLAB (Option 1)). Additional benefit is: 
It can extract the values of Field Output variables at all nodes or integration points at any prescribed load step.
This techniques is therefore very useful in extracting step-deformation data in simulations with multi-step loading.

Python scripts are kept in a separate subsubdirectory `python_scripts_for_odb_read`, although it is recommended to transfer these 
.py scripts to the C:\temp folder in your desktop so that you can use the default MATLAB scripts. 

# SECTION 5: POST PROCESSING

Examples of post processing used for my thesis and paper.

### Included Files
`ComputeYoungsModulus.m`

`GetResults.m`

### Description
`ComputeYoungsModulus.m` is a stripped down version of `GetResults.m`. The latter is a
way of getting computing all the information we used in computing bending and torsion
and tension of the same network and then computing. It offers some insight on to
everything I did but the prior is a good entryway in how to compute certain aspects of
the FE software.

In general, you'll have to customize the postprocessing to match your specific objective.

### Notes
Good luck!













