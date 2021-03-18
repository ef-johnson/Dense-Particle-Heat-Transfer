# This code reads in the pos file from LIGGGHTS, then finds the F matrix from the correlation (instead of reading it in from the
# output of the MC simulation, the normal way)

# Version 1.1 - vectorizing the PP_dist array calculation
# Version 1.2 - putting things into separate, external functions. Check into timing of spline2d, and how to pass interps into functions? 


using DelimitedFiles
using CSV
using Formatting
using Printf
using Dierckx

include("PP_ht.jl")
include("PW_ht.jl")
include("create_BCs.jl")
include("PFP_conduction.jl") 
include("PFW_conduction.jl") 
include("bilin_interp.jl")


function go()


# BASIC PARAMETERS
rebuild_BCs = true
rebuild_PFP_cond = true


LTS_timestep = 0.000002 # LIGGGHTS time step size, in seconds, from liggghts .in file

# Start and End time steps
LTS_stamp_delta =  2000 # nuumber of DEM time steps per thermal time step
LTS_stamp_start = 80000 # start at this DEM time step
LTS_stamp_end =  1000000 # end at this DEM time step

thermal_ts = LTS_timestep * LTS_stamp_delta # Thermal time step, seconds


println("Thermal time step size: $thermal_ts sec")
println("DEM time steps per thermal time step: $LTS_stamp_delta")


n_threads = Threads.nthreads() # number of threads is set in the .bashrc file (in Ubuntu) by adding this command "export JULIA_NUM_THREADS=4", or by using the same command in a terminal.
println("Number of threads: $n_threads")




# PARTICLE PROPERTIES

# This code currently works for a volume fraction that is precalculated and can be divided into the bulk and near wall regions. VlF_bulk used in PP and PW radiation, VlF_near_wall used only in PP radiation model, and VlF_bulk is used as a parameter in the PFP conduction model. These volume fractions must be precalculated, and the near-wall region is taken as within 5 radii of the wall (see PP_ht.jl)
VlF_bulk = 0.60
VlF_near_wall = 0.60

k_p = 2.0 # Thermal conductivity of particle material, W/mK
k_w = 16.6 # thermal conductivity of wall, W/mK


# The values for emissivity to be used in the RDF radiation lookup tables 
e_p = 0.86 
e_w = 0.60

# The emissivity used in the acutal PP radiation heat transfer equation. 
e_p_eqn = 0.86 

# Young's mod of particles
Y_i_DEM = 1e8 # Youngs modulus used in DEM simulation (typically much less than actual, to allow for a larger time step in DEM)
Y_i_real = 52e9 # Youngs modulus of the actual particles

# Young's mod of wall
Y_wall_DEM = 1e8 # Youngs mod, used in DEM (stainless steel)
Y_wall_real = 180e9 # Stainless steel

nu_p = 0.25
nu_wall = 0.30

# Factor to adjust for the different contact area in the DEM vs. actual cases. Interpreted from Zhou, Yu & Zulli, 2010, "A new computational method..." 
c_factor_PP = (Y_i_DEM/Y_i_real)^(1/5) 
c_factor_PW = ( ( ((1-nu_p^2)/Y_i_real) + ((1-nu_wall^2)/Y_wall_real) ) / ( ((1-nu_p^2)/Y_i_DEM) + ((1-nu_wall^2)/Y_wall_DEM) ) )^(1/5) 
println("c_factor_PP: $c_factor_PP")
println("c_factor_PW: $c_factor_PW")

# Cutoff distances: ignore any heat transfer beyond this number of radii
PP_cutoff = 7 # paticle-particle cutoff, radii
PW_cutoff = 7 # paticle-wall cutoff, radii

# For flowing tube problems with a periodic boundary condition, re-assign the particle temperature to T_inlet for any particle at a location above 
T_inlet = 293 # Assign particles this temp at inlet
y_reset_temp = 0.19 # Reset temp of any particle above this to the inlet temp (for problems with periodic BC which loops particles from bottom back to the top)


radius = 0.0005 #  particle radius in meters
sigma = 5.670373e-8 # Stefan-Boltzmann constant, W/m2K4
area = 4*pi*radius^2 # surface area of a single sphere
density = 2200 # density, kg/m^3
Cp = 1130 # specific heat, J/kgK
mass = (4/3*pi*radius^3)*density # particle mass


# DELETE OLD DATA
rm("../post", recursive=true)
mkpath("../post/xyzT")
mkpath("../post/q_PW_tot")
mkpath("../post/q_PW_cells")



# READ RDF TABLES FOR RADIATION
# Tables have columns of volume fraction of .25, .35, .45, .55 and .64, and rows for different PP and PW distances, in units of radii
 
# PP RDF
PP_table_name = string("../PP_RDF_tables/PP_RDF_table_ep_", format(e_p, precision=2), ".txt")
PP_RDF_table_DF = CSV.read(PP_table_name, header=false,delim=' ')

# Set up PP arrays for interpolation
PP_RDF_table_arr = convert(Matrix, PP_RDF_table_DF) # Convert the table from a DataFrame (how CSV.read reads in data) to a normal array 
PP_RDF_dist =  PP_RDF_table_arr[2:end,1]
PP_RDF_VlF = PP_RDF_table_arr[1,2:end]
PP_RDF_values = PP_RDF_table_arr[2:end,2:end] 


# If cutoff is set higher than the table, reset the cutoff to the end of the table (avoids extrapolation)
if PP_cutoff > maximum(PP_RDF_dist)
	PP_cutoff = maximum(PP_RDF_dist)
end
println("PP cutoff: $PP_cutoff")


# PW RDF
PW_table_name = string("../PW_RDF_tables/PW_RDF_table_ep_", format(e_p,precision=2), "_ew_",format(e_w,precision=2), ".txt")
PW_RDF_table_DF = CSV.read(PW_table_name, header=false, delim=' ')

# Set up PW arrays for interpolation
PW_RDF_table_arr = convert(Matrix, PW_RDF_table_DF) # Convert the table from a DataFrame (how CSV.read reads in data) to a normal array 
PW_RDF_dist =  PW_RDF_table_arr[2:end,1]
PW_RDF_VlF = PW_RDF_table_arr[1,2:end]
PW_RDF_values = PW_RDF_table_arr[2:end,2:end] 


if PW_cutoff > maximum(PW_RDF_dist)
	PW_cutoff = maximum(PW_RDF_dist) # Keep this high for now, maybe want to reduce it later
end
println("PW cutoff: $PW_cutoff")






# PFP and PFW HTC TABLES

# Heat transfer coefficients are precalculated by numerical integration in the PFP_conduction function. 
# The distances and temps to calculate are defined in the function. 
# The function writes a matrix of htc values at various temperatures and PP distances, which are written into PFP_dist.txt and PFF_Temps.txt.
 
# For PFP, the distance only needs to go up to 3*radius per recommendations by XXX researchers, but for PFW, a ghost particle is used, so if you want 
# the distances to be calculated up to 1 radius between tip of particle and the wall, that is the same as 2 radii between the particle and its ghost particle
# so a maximum of 4 radii center-to-center is needed. 

if rebuild_PFP_cond == true
	PFP_conduction(VlF_bulk, k_p, radius)
	PFW_conduction(VlF_bulk, k_p, radius)
end 

# Read data output by function PFP_conduction
PFP_all = CSV.read("../PFP_PFW_conduction/PFP_htc_values.txt", header=false, delim=' ') # Reads in text files as a DataFrames
PFP_all = convert(Matrix,PFP_all) # Convert from DataFrame to normal matrix
PFP_dist = PFP_all[2:end,1]
PFP_Temps = PFP_all[1,2:end]
PFP_htc_values = PFP_all[2:end,2:end]

# Read data output by function function PFW_conduction
PFW_all = CSV.read("../PFP_PFW_conduction/PFW_htc_values.txt", header=false, delim=' ') # Reads in text files as a DataFrames
PFW_all = convert(Matrix,PFW_all) # Convert from DataFrame to normal matrix
PFW_dist = PFW_all[2:end,1]
PFW_Temps = PFW_all[1,2:end]
PFW_htc_values = PFW_all[2:end,2:end]



# CREATE BOUNDARY CONDITIONS

# Read in wall and boundary condition file. If rebuild_BCs is true, this function will run to create the BC file, bdry_data
# Boundary conditions must be set in the create_BC function file.
# Output of function create_BCs() is an array with columns for each wall element (aka cell): Cell Center X, Cell Center Y, Cell Center Z, a, b, c, d, Temp, BC_type, Area
# Cell centers X, Y, Z are the centroid locations, the coefficients a to d are for the plane passing of each cell of form ax + by + cz + d = 0, which is used used to find
# the perpendicular PW distance). Temp is temperature. BC_type is the type of the boundary condition where 1 is a specified temperature, and 2 is adiabatic (aka heat transfer 
# from particle to this wall element is 0)

if rebuild_BCs == true
	create_BCs()
end 

# Read in the BC data
bdry_data = CSV.read("../mesh/bdry_data.txt", header=false,delim=' ')
bdry_data = bdry_data[:,1:10]
bdry_data = convert(Matrix, bdry_data)

centers = bdry_data[:, 1:3]
plane_coeffs = bdry_data[:,4:7]
T_wall = bdry_data[:,8]
BC_type = bdry_data[:,9]
element_area = bdry_data[:,10]
n_cells = size(bdry_data,1)
println("Number of wall elements: $n_cells")






# --------------------------TIME STEP LOOP-------------------------

ts_count = 0

for LTS_stamp = LTS_stamp_start : LTS_stamp_delta : LTS_stamp_end
t_beg = time_ns()

ts_count = ts_count + 1
real_time = (ts_count-1) * thermal_ts # Time of heat transfer simulated



# READ CURRENT POSITION FILE
# position file should have be names pos_[timestep].txt and have columns x y z, space delimited.
pos_name = string("../DEM/post_positions/pos_", format(LTS_stamp),".txt")
pos = CSV.read(pos_name, header=false,delim=' ')
#pos = CSV.read(pos_name, header=false, datarow=10, delim=' ') # Use this if using the default dump LIGGGHTS file format (skips the 9 lines of header)
pos = pos[:, 1:3] # sometimes an extra column of spaces is added from LIGGGHTS, so remove this (and any other data that might have been output, like velocities)
# println("pos_DF: ", typeof(pos)) # Can check the type with these lines. If it shows a data type with "missing", then code may still work, but very slowly.  
# display(pos)
pos = convert(Matrix, pos)
xPos = pos[:,1]
yPos = pos[:,2]
zPos = pos[:,3]

n_particles = size(pos,1) # number of particles in simulation - should be constant throughout simulation




# READ AND INITIALIZE TEMPS

T_curr = zeros(n_particles,1)

# Initialize Temperatures if it is the first time step
if ts_count == 1

	T_curr = ones(n_particles,1)*T_inlet
	println("Number of particles: $n_particles")	
	
# For all other time steps, read the xyzT file
else
	temp_name_prev = string("../post/xyzT/xyzT_", format(ts_count-1),".txt")
	T_curr = (readdlm(temp_name_prev, ' ', Float64, '\n'))
	T_curr = T_curr[:,4] # T_curr are the particle temps at the last time step 
	
end

T_next = zeros(n_particles,1) # temperature that will be solved for at this time step


# PP HEAT TRANSFER 

# Uses Threads.@threads to paralellize the for loops. Should be used inside the function for best speedup. 

t_PP_beg = time_ns()

q_PP = PP_ht(n_particles, T_curr, xPos, yPos, zPos, radius, PP_cutoff, k_p, c_factor_PP, VlF_near_wall, VlF_bulk, PP_RDF_dist, PP_RDF_VlF, PP_RDF_values, e_p_eqn, area, sigma, PFP_dist, PFP_Temps, PFP_htc_values)


t_PP_end = time_ns()
#display(q_PP)



# PW HEAT TRANSFER 

# Sum of heat transferred to each wall element from each ht mode
q_cond_wall = zeros(n_cells,1)
q_conv_wall = zeros(n_cells,1) 
q_rad_wall = zeros(n_cells,1)

# Function PW_ht() outputs the q_PW[i] array, and it also changes q_cond_wall, q_conv_wall, and q_rad_wall arrays from zeros to have the calculated values. 
q_PW = PW_ht(n_particles, VlF_bulk, radius, PW_cutoff, xPos, yPos, zPos, centers, n_cells, BC_type, plane_coeffs, c_factor_PW, T_wall, T_curr, k_p, k_w, PFW_dist, PFW_Temps, PFW_htc_values, PW_RDF_dist, PW_RDF_VlF, PW_RDF_values, e_p_eqn, area, sigma, q_cond_wall, q_conv_wall, q_rad_wall)

t_PW_end = time_ns()


 
# Total HT from particles to all walls
q_PW_tot = sum(q_PW) # Total heat transfer from all particles to walls
q_PW_cond_tot = sum(q_cond_wall) # Total conduction ht to walls
q_PW_conv_tot = sum(q_conv_wall) # Total convection ht to walls
q_PW_rad_tot  = sum(q_rad_wall) # Total radiation ht to walls

q_all_wall = q_cond_wall + q_conv_wall + q_rad_wall # Total ht to each wall element


# RECALCULATE PARTICLE TEMPS
for i = 1:n_particles

	# Calculate new particle temps for particles below y_reset_temp
	if yPos[i] < y_reset_temp
		T_next[i] = T_curr[i] +  (q_PP[i] + q_PW[i]) * thermal_ts / (mass * Cp)
	else
		T_next[i] = T_inlet
	end
				
end


# OUTPUTS FOR POST PROCESSING AND VISUALIZATION

# xyzT file - prints position and temperature of each particle at every time step
all_out_name = string("../post/xyzT/xyzT_", format(ts_count),".txt")

open(all_out_name, "w") do io
	for i = 1:n_particles
		@printf(io, "%.8f %.8f %.8f %.6f\n", xPos[i], yPos[i], zPos[i], T_next[i])
	end
end



# q_PW TOTAL FILE - columns: Time step number, total PW ht (from all modes), total PW cond ht, total PW conv ht, total PW rad ht 
q_PW_tot_name = string("../post/q_PW_tot/q_PW_tot.txt")
open(q_PW_tot_name, "a") do io
		@printf(io, "%.0f %.8f %.8f %.8f %.8f\n", ts_count, q_PW_tot, q_PW_cond_tot, q_PW_conv_tot, q_PW_rad_tot)
end

# q_PW_cells FILE: ht to each wall element, columns: centroid x, y, z, total ht to element, total cond ht, total conv ht, total rad ht. 
q_PW_name = string("../post/q_PW_cells/q_PW_cells_", format(ts_count),".txt")
open(q_PW_name, "w") do io
	for i = 1:n_cells
		@printf(io, "%.8e %.8e %.8e %.8e %.8e %.8e %.8e \n", centers[i,1], centers[i,2], centers[i,3], q_all_wall[i], q_cond_wall[i], q_conv_wall[i], q_rad_wall[i])
	end
end


# TIMING
t_end = time_ns()

t_init = (t_PP_beg - t_beg) / 1e9
t_post = (t_end - t_PW_end) / 1e9
t_loop = (t_end - t_beg) / 1e9
t_PP = (t_PP_end - t_PP_beg) / 1e9
t_PW = (t_PW_end - t_PP_end) / 1e9

t_init = format(t_init, precision = 2)
t_post = format(t_post, precision = 2)
t_loop = format(t_loop, precision = 2)
t_PP = format(t_PP, precision = 2)
t_PW = format(t_PW, precision = 2)

q_PW_tot = format(q_PW_tot, precision=3)

real_time = format(real_time, precision=4)

println("Time step: $ts_count, DEM ts: $LTS_stamp, real time: $real_time, q_PW total: $q_PW_tot \n	Loop timing: Total = $t_loop, Init = $t_init, PP =  $t_PP, PW = $t_PW, post = $t_post ")


end # end time step loop



end; # end of function

