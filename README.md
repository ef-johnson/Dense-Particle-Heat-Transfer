# Dense-Particle-Heat-Transfer
This code calculates heat transfer through groups of small, sand-like particles, in dense granular flows. It was developed for particle heat exchangers and solid particle solar receivers, but it may also be useful for research on similar topics, such as lime kilns, laser sintering of powder beds, or metallurgical processes. 

DPHT is a heat transfer code for large groups of small spherical particles, which are either static or in dense granular flows. It is based on the Discrete Element Method (DEM). The particle collision mechanics are first simulated with DEM, and the particle xyz positions are written to text files over time. After the DEM simulation is finished, this DPHT code reads each of these xyz position files and calculates the heat transfer between the particles, as well as between particles and the walls. This can be considered a one-way coupled to DEM. 

DPHT was developed for modeling solid particle solar receivers and heat exchangers in the field of concentrating solar power, but it is applicable for many situations having packed or moving particle beds. The heat transfer calculations are based on work from previously published research, which is explained in the thesis:
Advances in Modeling High Temperature Particle Flows in the Field of Concentrating Solar Power, Evan F. Johnson, Middle East Technical University, 2021. 
Chapter 5 of the thesis contains the background, assumptions, model details, and an example. The code for the same example is included in this GitHub repository.  

It is recommended to work in Ubuntu or a similar Linux-based operating system, so that other related software works easily (especially LIGGGHTS for DEM, and ParaView for visualization).  

Prerequisites:

A) The DPHT code is written in Julia, which is extremely fast due to its run-time compilation, but coding it is much simpler than working in the other compiled languages like C/C++/Fortran. If you know Matlab or similar languages, Julia is easy to adopt. Julia can be downloaded for free (https://julialang.org/). 

B) Before DPHT can be run, it needs a DEM code to solve for the particle positions over time. The open source code LIGGGHTS is recommended (https://www.cfdem.com/media/DEIt is presented "as-is", and no guarantees are made to its accuracy. Extensive notes are written in the comments to help others interpret the codM/docu/Manual.html). Within the "DEM" directory of the DPHT file structure, the LIGGGHTS input file (the "in" file) is stored. When running the DEM simulation, the particle center positions are output in two places: A) The folder "post_positions", where the file format is a simple text file (these files will be read by DPHT), and B)The "post_vtk" folder, which is slightly easier to view in ParaView. 

C) ParaView is needed for visualizing the particle positions, temperatures, etc. It is also used to find the wall element centers and the relations between triangular mesh elements. 

D) Some type of CAD drawing software is needed for constructing the geometry that contains the particles, such as a tube, a heat exchanger surface, etc. 

E) A meshing program is needed which can make 2D surface meshes with triangular elements, such as Gmsh.  


To run a simulation, the steps are given below. Much more detail is given in the thesis. 

1) Draw the geometry of the walls containing the particles using a CAD program (SolidWorks, Onshape, etc).
2) Mesh the geometry with triangular elements (using a program such as Gmsh) and save the mesh in STL file format in the "mesh" folder.
3) Read the STL file in to ParaView, where three text files must be created and saved in the "mesh" folder: A)  vertices.txt, which contains the coordinates of the points of the triangles, B) centers.txt, which contains the center coordinates of the triangles, and C) relations.txt, which gives the relationship between the element IDs and the vertex IDs.
5) The DEM simulation is run using LIGGGHTS, from within the "DEM" folder. This requires opening a terminal and issuing a command such as "mpirun -np 20 ~/LIGGGHTS/LIGGGHTS-PUBLIC/src/lmp_auto < in.fill_lattice" (assuming LIGGGHTS has been downloaded and the executable "lmp_auto" has been created. (Note that one change should be made before compiling LIGGGHTS: the output ("dump") files should have their header lines removed for simplest use - see thesis.)
6) The parameters for the heat transfer simulation can be set in the "DPHT_code" folder. In the tube example given, the file "DPHT_tube.jl" contains the main function with the thermal properties, time step information, initial particle temperatures, etc. The wall temperatures are set in the "create_BCs.jl" file. Once the DEM simulation is finished, open a terminal window in the "DPHT_code" folder, start Julia, and "include" the various files with the command: include("DPHT_tube.jl"). Then the simulation is started by calling the main function of DPHT_tube.jl with the command: "go()".
7) DPHT outputs are saved in the "post" folder. The particle xyz positions and temperatures can be visualized with ParaView using the "xyzT[timestep].txt" files - one file is written out for each time step. The "q_PW_tot.txt" file shows the total particle-wall heat transfer at each time step, which is useful to know for a heat exchanger or for judging when a system has reached steady state. The "q_PW_cells.txt" files give the particle-wall heat transfer for each triangular wall element, at each time step.

The code is my first time writing in Julia, so there are undoubtedly ways to improve it. Extensive notes are written in the comments to help others interpret the code. It is presented "as-is", and no guarantees are made to its accuracy.

This code was developed as part of the mentioned thesis, but it is now shared so others can use it and build upon it. If you use this code or find the thesis useful for your own work, please cite it. Thanks!








