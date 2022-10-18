# Readme for Network Generation

*This repositry has methods to generate fiber networks. A description is below.*

## Useful Resources to Read
* Peter Grimmer's Thesis
* Grimmer, P. and Notbohm, J. Displacement propagation in fibrous networks due to local contraction. J Biomech Eng, 140(4):041011, 2018.
* Stephen Tyznik's Thesis
* Lindstrom et al, 2010, Phys Rev E, 82, 5, 051905; Lindstrom et al, 2013, Soft Matter, 9, 30, 7302

## Broad Overview

There are different subdirectories, each containing a different method to generate a fiber network.

* simulated_annealing_3D: This method makes 3D networks matching the structure of collagen gels, and is based on the manuscripts by Lindstrom et al.

* simulated_annealing_2D: This is a 2D version of the 3D simulated annealing.

* lattice_2D: This directory contains methods to make 2D networks from nodes originally arranged in a lattice. This method creates networks that don't match collagen gels perfectly, but these methods are commonly used in other manuscripts.


## Specific Details About the Different Methods

### Simulated Annealing (3D)

This method was first proposed in two manuscripts by Lindstrom et al (Lindstrom et al, 2010, Phys Rev E, 82, 5, 051905; 
Lindstrom et al, 2013, Soft Matter, 9, 30, 7302). Our version of the code was first implemented by Peter Grimmer and 
then later improved by Stephen Tyznik.

Use `Fiber_Network_Generation_Main_Script.m`

There are other subfunctions; read comments in those subfunctions along with theses by Peter Grimmer and 
Stephen Tyznik for more details.

There exists one artifact in these simulated networks. This algorithm has a tendency to deposit more nodes near the edges of the x--y planes comprising the 3D network domain. To achieve an uniform nodal density 
throughout the network domain, a few nodes and fibers sufficiently close to these edges should be cropped out. The function file `auto_crop_tool.m` and the subfunction file `check_raw_cropped.m` 
provided within the folder fiber_network_model/1-NetworkGeneration/simulated_annnealing_3D/artifact_remover_3D_network does this job. In general, one should remove fibers and nodes lying within 
a distance of 15% of the respective domain dimension from the edges (as viewed in plan) by inputting p = 0.7 in the function `auto_crop_tool.m`.
 
### Simulated Annnealing (2D)

This subdirectory modified the 3D simulated annealing code to generate a 2D network using the simulated annealing method.
Code in this subdirecotry was originally written by Mainak Sarkar. Mainak modified the simulated annealing-based branch 
optimization steps in all the function files and main script from 3D to 2D. 

##### Main script

`Fiber_network_main_generation_2D.m`

##### Supporting function files

`Initial_Network_Generation_2D.m`

`Network_Optimization_2D.m`

`Branch_Optimization_2D.m`

##### Description
Running the script, `Fiber_network_main_generation_2D.m`, is all that needs to be
done to create the 2D simulated annearling fiber network geometry to be used in Finite Element applications. This
script will call the remaining functions in order and then save the model in a folder with its name.

In `Fiber_network_main_generation_2D.m` you are able to change network variables
such as dimensions, nodal density, RNG seed, and fiber length. Additionally, you are
able to set the save name for the network. Intelligently name your files based on the
variables you assign it so it is easy to find correct networks for future use.
`Initial_Network_Generation_2D.m` creates a network that satisfies the connectivity
distribution defined as well as the nodal density. The fiber lengths are optimized
through an iterative process in `Network_Optimization_2D.m`. Similarly branching angles are
optimized in the final function, `Branch_Optimization_2D.m`.

##### Notes
When altering properties such as density or fiber length. Double check distribution of
important variables as it is not always clear how they are coupled. For example, when
increasing fiber length, artificial fiber alignment was generated in my experiments.



### Lattice (2D)

There are various options for different methods that start with a lattice of nodes:
* Triangular lattice-based
* Hexagonal polygon-based 
* Voronoi tesselation-based 

All generate 2D networks of fibers.

Code in this subdirectory was originally written by Mainak Sarkar.

##### Useful references

*Ronceray P., Broedersz C. P, Lenz M., Fiber networks amplify active stress. P Natl Acad Sci USA, 2016, 113, 2827.


##### Main script to generate the fiber network

`lattice_main_network.m`

##### Function files and brief info

1.	`add_midnode_in_fiber.m` - Applicable to other networks 
2.	`arbit_node_neighbors.m` - Specific to lattice networks 
3.	`dangling_fiber_remover.m` - Applicable to other networks 
4.	`delaunay_tri_mesh_gen.m` - Specific to lattice networks 
5.	`edge_element_generator.m` - Specific to lattice networks 
6.	`hexa_mesh_gen.m` - Specific to lattice networks 
7.	`nodal_disorder.m` - Specific to lattice networks 
8.	`random_fiber_eliminator.m` and `random_fiber_eliminator_adv.m` - Applicable to other networks 
9.	`renumbering_el_nodes.m` - Applicable to other networks 
10.	`find_node_numbers_2Dlattice.m` - Applicable to 2D lattice networks only
11. `ghost_node_remover.m` - It is also a very general function file, applicable to any network
12. `find_avg_node_connectivity.m` - Useful for any network.
13. `draw_edge_fibers.m` - It is best suitable for polygonal (Voronoi) lattice-based networks to make robust loading edge
14. `rhonodal_finder.m` - For 2D lattice network only
15. `distance_finder.m` - Applicable to other networks 

##### Description of generated network:

Code addresses two methods to generate fiber network similar to what is described in the manuscript by Ronceray et al.
There is a switch to choose method 1 or method 2 at the beginning. (Note: In the following description, even though 
I am saying triangular lattice, it actually means there is also a mid-node at each fiber.) 

Method 1: The fundamental lattice is triangular. The constraint is that each node is connected to exactly six fibers. 
Under this method 1, I have a switch to derive a hexagonal network from this, with the specific constraint that 
each node is connected to exactly 3 fibers. Method 1 does not form a DT or a Voronoi network. 

Method 2: The fundamental lattice is again triangular. But the network is obtained through Delaunay triangulation 
applied on the random nodes in the 2D domain. So, here, I do not have a constraint of 6 fibers to be connected to 
one node. Delaunay triangulation maximizes the minimum angle of all the angles of the triangles in the triangulation, 
so I have introduced this method is my code. This should be the most fundamental lattice network (and compact too) 
as far as Ronceray’s paper is concerned.  Further, the theory says, 

>The Delaunay triangulation of a discrete point set corresponds to the dual graph of the Voronoi diagram of that point set. The circumcenters of Delaunay triangles 
are the vertices of the Voronoi diagram. In our case, the Voronoi vertices are connected via edges, that can be 
derived from adjacency-relationships of the Delaunay triangles: If two triangles share an edge in the Delaunay 
triangulation, their circumcenters are to be connected with an edge in the Voronoi tesselation.

So, I have applied this Voronoi tesselation to obtain the corresponding Voronoi network from the Delaunay triangulated network obtained in method 2. 
So, in brief, I have a switch in method 2 in the main script, to choose either only a Delaunay triangulated network or a (DT + derived Voronoi) polygonal network.   

Other operations on the network, including removing random fibers based on probability of existence of each fiber in the network,
removing dangling fibers, adding mid-nodes to each fibers, and other function files work on both methods. 

##### Comprehensive input panel

There are many binary switch-based choices to customize the flow of code. The input panel describes the entries.

##### A brief overview of the steps in MATLAB code 

1.	The 2D lattice / Delaunay / Voronoi tessellation-based fiber network is generated in these codes. Then an Abaqus input file (INP) of the network is created.  
2.	Abaqus input file (INP) is run within the MATLAB through DOS to generate the output database file in ODB format.  
3.	The required field output variables are extracted from this ODB file through Python scripts, run via MATLAB. Output variables are stored in a MAT file. 

The important significance of steps 2 and 3 above are: 

a)	Fully controllable from MATLAB environment. 

b)	We do not need to open Abaqus GUI to perform our simulation and extract output results. 

c)	In future simulations, if we need to run a loop of INP files in Abaqus, we can add a loop. 



 