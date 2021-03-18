# Function to set boundary conditions.
# Requires inputs of the nodes, cell centers, and cell-node relation info, as output by paraview. These must be named nodes.txt, cell_ctr.txt, and relations.txt.
# Outputs a text file with the necessary boundary condition info, bdry_data.txt


function create_BCs()

println("\nExecuting create_BC_array function. \nCreating boundary conditions and exporting file to mesh folder.")


# READ IN THE TEXT FILES OUTPUT BY PARAVIEW

# Nodes will be referred to here as their row number, so the node ID goes from 1 to n_nodes
#nodes_name = "../mesh/nodes.txt"
#nodes = CSV.read(nodes_name, header=true,delim=',')
#nodes = convert(Matrix, nodes) # convert to an array
#n_nodes = size(nodes,1)
#display(nodes)

vertices_name = "../mesh/vertices.txt"
vertices = CSV.read(vertices_name, header=true,delim=',') # This should have columns of: Points_0,Points_1,Points_2,Points_Magnitude,Point ID
vertices = convert(Matrix, vertices) # convert to an array
n_vertices = size(vertices,1)

vertices = sortslices(vertices,dims=1,by=x->x[5],rev=false) # sorts on the Point ID (in case it was output from Paraview sorted in a different order)
vertices = vertices[:,1:3] # Just take the xyz positions, discard the Points_Magnitude and Point ID columns


# Relations identify the thee vertices (points of the triangle) associated with each element center.
# By default, after "Export Scene" from paraview, there are 6 columns: Cell Type,Point Index 0,Point Index 1,Point Index 2,STLSolidLabeling,Cell ID
# First, sort based on the Cell ID, in case it was output after sorting by something else (but Cell ID is the default thing to sort by, so usually it's already sorted)
# Here, the Type, STLSolidLabeling and Cell ID we will discard. The Cell IDs range from 0 up to n_elements-1. But once sorted, the row number becomes the de-facto ID, so it's then 1 to n_elements. 
# The vertices are numbered 0 to n_vertices-1 by default, but we want them to range from 1 to n_vertices, so they match the row number used in "vertices" array - so add 1 to all the ids in the relations array.
relations_name = "../mesh/relations.txt"
relations = CSV.read(relations_name, header=true,delim=',')
relations = convert(Matrix, relations) 
relations = sortslices(relations,dims=1,by=x->x[6],rev=false) # sorts on the Cell ID (in case it was output from Paraview sorted in a different order)

relations = relations[:,2:4] # Just take the numbers of the three nodes
relations = relations .+ 1 # Have to add 1 to each of the node IDs, so they start at 1 instead of 0.
n_elements = size(relations,1)


# Element centers. Should have columns: Points_0,Points_1,Points_2,Points_Magnitude,STLSolidLabeling,Point ID
centers_name = "../mesh/centers.txt"
centers = CSV.read(centers_name, header=true,delim=',')
centers = convert(Matrix, centers) # turn it into an array
centers = sortslices(centers,dims=1,by=x->x[6],rev=false) # sorts on the Point ID (in case it was output from Paraview sorted in a different order)
centers = centers[:,1:3] # Just take the xyz positions, discard the Points_Magnitude, STLSolidLabeling, and Point ID columns
n_centers = size(centers,1)


println("$n_elements elements, and $n_vertices vertices\n")

if n_elements != n_centers
println("n_elements doesn't equal n_centers! Problem.")
end


bdry_data = []

# Iterate through all elements
for i = 1:n_elements
	
	# Get the xyz position of the cell center
	xCell = centers[i,1]
	yCell = centers[i,2]
	zCell = centers[i,3]
	
	
	# ASSIGN BOUNDARY CONDITIONS
		
	# BC_type: 1 for first kind (Temp specified), 2 for adiabatic
	
	# First set all to type 2 BC, and any temp (doesn't matter), then reassign certain elements with temp BC.
	BC_type = 2
	Temp = 99
	
	# For hot and cold walls, reassign them to the hot and cold temps, and BC type 1.
	if yCell > 0 && yCell < 0.14 # Set temp for heated part of the wall
		Temp = 1073 - yCell*1000 # Kelvin. It's 1073 at y = 0 and 923 at y=.015
		BC_type = 1
	end
	
	
	# Find the coefficients of the plane formed by each triangular cell.
	
	# Find the three node id's for this cell
	n1_id = relations[i, 1]
	n2_id = relations[i, 2]
	n3_id = relations[i, 3]
	
	# Look up the xyz locations of these nodes
	pos1 = vertices[n1_id, :]
	pos2 = vertices[n2_id, :]
	pos3 = vertices[n3_id, :]
	

	# Calculate the plane coefficients, in form ax + by + cz + d = 0.Cross prod of two vectors gives an equation for the line normal to both, 
	# aka normal to the plane formed by both. A plane has the equation ax + by + cz + d = 0, where a,b,c is a vector normal to the plane, 
	# so after taking the cross product we already know three of the coefficients. Plug in any point to find d.
	
	AB = pos2 - pos1 # Vector from node 1 to 2
	AC = pos3 - pos1 # Vector from node 1 to 3
	
	# Cross product to find normal vector (a,b,c)
	a = AB[2]*AC[3] - AB[3]*AC[2]
	b = AB[3]*AC[1] - AB[1]*AC[3]
	c = AB[1]*AC[2] - AB[2]*AC[1]
	
	# Plug in one point (doesn't matter which one) to find d
	d = -(a*pos1[1] + b*pos1[2] + c*pos1[3])


	# Find the surface area of the element. Area = (1/2)*magnitude(cross product)
	# We already have the cross product vector calculated as a,b,c
	Area = (1/2) * sqrt(a^2 + b^2 + c^2)
	

		
	# Accumulate data in an array with columns: Cell Center X, Cell Center Y, Cell Center Z, a, b, c, d, Temp, Area 	
	append = [centers[i,1] centers[i,2] centers[i,3] a b c d Temp BC_type Area]
	if i==1
		bdry_data = append
	else	
		bdry_data = vcat(bdry_data, append)
	end
	

end
	
	# Write to file bdry_data.txt
	writedlm("../mesh/bdry_data.txt", bdry_data, ' ')

end




























