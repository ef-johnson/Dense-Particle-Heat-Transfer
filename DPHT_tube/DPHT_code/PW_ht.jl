# Calculates the net heat transfer rate to each particle from all wall elements, q_PW. Positive values indicate heat transferred TO a particle, negative is heat FROM a particle.
# Returns q_PW to main code, and changes the vectors q_PW_cond, q_PW_conv, q_PW_rad to be used in the main code.  

function PW_ht(n_particles, VlF_bulk, radius, PW_cutoff, xPos, yPos, zPos, centers, n_cells, BC_type, plane_coeffs, c_factor_PW, T_wall, T_curr, k_p, k_w, PFW_dist, PFW_Temps, PFW_htc_values, PW_RDF_dist, PW_RDF_VlF, PW_RDF_values, e_p_eqn, area, sigma, q_cond_wall, q_conv_wall, q_rad_wall)


q_PW = zeros(n_particles,1)

Threads.@threads for i = 1:n_particles

	
#Initialize ht from particle i to wall as zero
q_PW_cond = 0
q_PW_conv = 0
q_PW_rad = 0
	
											
# For PW heat transfer, use VlF_bulk
VlF = VlF_bulk

# Check if particle is close to a wall element: write a FOR loop so that PW heat transfer isn't evaluated
# if the particle is a long ways from the all. Use PW_cutoff for this. 
	
# For example, for parallel plates where the heated wall is at z=0, use:
# "if zPos[i]/radius < PW_cutoff" 
# so the PW heat transfer is only evaluated for particles that are close to the wall. 
# If not, PW heat transfer will be left as 0. Just to save computation time. 
	

	
	# Find the closest wall element center
	xPos_i = xPos[i]
	yPos_i = yPos[i]
	zPos_i = zPos[i]
		
	curr_min = Inf # current minimum distance
	curr_min_id = -99 # current minimum wall element id
	
	for j = 1:n_cells
	
		P_cell_center_dist = ( (xPos_i - centers[j,1])^2 + (yPos_i - centers[j,2])^2 + (zPos_i - centers[j,3])^2 )^(0.5)
		
		if P_cell_center_dist < curr_min
			curr_min = P_cell_center_dist
			curr_min_id = j
		end
		
	end
			
	wall_id_closest = curr_min_id	
	
	# Only compute PW heat transfer if the closest cell has a type 1 BC
	if BC_type[wall_id_closest] == 1
		
		# Use the pre-calculated plane coefficients for each wall element, in the form of ax + by + cz + d = 0
		coeff_a = plane_coeffs[wall_id_closest,1]
		coeff_b = plane_coeffs[wall_id_closest,2]
		coeff_c = plane_coeffs[wall_id_closest,3]
		coeff_d = plane_coeffs[wall_id_closest,4]	
						
		# Find the distance from particle center to plane of the wall element which has the closest center	
		# Formula for point to plane distance: 
		# https://mathworld.wolfram.com/Point-PlaneDistance.html						
		PW_dist = abs( (coeff_a*xPos[i] + coeff_b*yPos[i] + coeff_c*zPos[i] + coeff_d) / sqrt(coeff_a^2 + coeff_b^2 + coeff_c^2) )  
		PW_dist_norm = PW_dist/radius
	
	
	
		# PARTICLE WALL CONDUCTION			
		if PW_dist_norm < 1 # Particle is overlapping the wall
			r_c = sqrt(radius^2 - PW_dist^2) # radius of contact
				
			q_PW_cond = 4 * r_c * c_factor_PW * (T_wall[wall_id_closest] - T_curr[i]) / (1/k_p + 1/k_w) 
			
			# Add this conduction to the total for the corresponding wall element
			q_cond_wall[wall_id_closest] = q_cond_wall[wall_id_closest] + q_PW_cond
		end
		
		
		
		
		# PARTICLE FLUID WALL CONDUCTION (aka PW Convection)
									
		if (PW_dist_norm < 1.5) # only proceed H/R < 0.5
						

#			PP_dist_ghost = 2*PW_dist
			T_eval =	(T_wall[wall_id_closest] + T_curr[i]) / 2 # Evaluate the PFW htc at the mean temp of the particle and the wall		
			
			# Correct for reduced Young's mod if particle is overlapping the wall			
			if PW_dist_norm < 1 # if overlapping
				PW_dist_real = sqrt( radius^2 - c_factor_PW^2 * (radius^2 - PW_dist^2) )
			else
				PW_dist_real = PW_dist # not overlaping
			end
			# bilin_interp(row values, col values, data, row_i, col_i, function name so that this bilin_interp is identified if there's an error)					
			htc_PFW = bilin_interp(PFW_dist, PFW_Temps, PFW_htc_values, PW_dist_real, T_eval, "PW convection") 
			q_PW_conv = (T_wall[wall_id_closest] - T_curr[i]) * htc_PFW
			
			# Add this convection to the total for the corresponding wall element
			q_conv_wall[wall_id_closest] = q_conv_wall[wall_id_closest] + q_PW_conv
	
		end # end if particle is within 2 radii of wall	
				

				



		# PARTICLE WALL RADIATION
			
		# If particle overlaps the wall, reset to 1, per DBA method
		if PW_dist_norm < 1.0
			PW_dist_norm = 1.0
		end
			
		# Find the proper PW RDF with a bilinear interpolatio nover distance and volume fraction							
		PW_RDF = bilin_interp(PW_RDF_dist, PW_RDF_VlF, PW_RDF_values, PW_dist_norm, VlF, "PW radiation") 

		# PW from particle i to closest wall element
		q_PW_rad =  PW_RDF * e_p_eqn * area * sigma * (T_wall[wall_id_closest]^4 - T_curr[i]^4) 

		# Add this radiation to the total for the corresponding wall element
		q_rad_wall[wall_id_closest] = q_rad_wall[wall_id_closest] + q_PW_rad	 # add to wall element rad total

	end # end if the closest wall element has BC type 1
		
#end # end if particle is within PW_cutoff of a wall element


q_PW_tot_i = q_PW_cond + q_PW_conv + q_PW_rad # May be faster for multithreading this way
	
q_PW[i] = q_PW[i] + q_PW_tot_i
	

end # iterate through i particles

return q_PW

end # end function

