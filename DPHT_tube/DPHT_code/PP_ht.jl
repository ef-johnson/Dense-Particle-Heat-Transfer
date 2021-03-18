# Calculates the net heat transfer rate to each particle, q_PP. Positive values indicate heat transferred TO a particle, negative is heat FROM a particle.
# Returns q_PP to main code. 

function PP_ht(n_particles, T_curr, xPos, yPos, zPos, radius, PP_cutoff, k_p, c_factor_PP, VlF_near_wall, VlF_bulk, PP_RDF_dist, PP_RDF_VlF, PP_RDF_values, e_p_eqn, area, sigma, PFP_dist, PFP_Temps, PFP_htc_values)


# For optimized distribution of work to each thread, build a vector "interleaved_ids", which alternates between particle id (i) values near the beginning and near the end.
# This is used in q_PP since low i's have a low workload compared to high i's, due to the for loop structure where j runs from 1 to (i-1). 
interleaved_ids = Vector{Int64}(undef, n_particles)

cnt_up = 1
cnt_down = n_particles
for l = 1:n_particles
	if l%2==1 # if i is odd
		interleaved_ids[l] = cnt_up
		cnt_up = cnt_up + 1		
	else
		interleaved_ids[l] = cnt_down
		cnt_down = cnt_down - 1
	end	
end

q_PP = zeros(n_particles,1)

# Start multithreaded for loop.
Threads.@threads for k = 1:n_particles

i = interleaved_ids[k] # reassign the i's so they're balanced

T_curr_i_4 = T_curr[i]^4 # Precalculate T^4 for particle i, so it doesn't have to recalculate it for every j iteration

	# Iterate through the other particles, for this i particle. Since ht from i to j is same as from j to i, only need to iterate from j = 1 to i. 
	# For j=i, this is ht to self, which is zero, so iteration is from j = 1:(i-1). Once ht is calculated for each ij combination, heat is added
	# to particle i and subtracted from particle j. Compared to the basic option of iterating through all j's, this is 2x faster, and certain to conserve energy.  

	for j = 1:(i-1)
		
		# Find PP distance, contact radius
		PP_dist = sqrt((xPos[i]-xPos[j])^2 + (yPos[i]-yPos[j])^2 + (zPos[i]-zPos[j])^2 )  # PP distance in meters. 		
		PP_dist_norm = PP_dist / radius # normalized by radius
		
		if PP_dist_norm < PP_cutoff # Only continue if it's below the cutoff distance
			
		
			# PP CONDUCTION				
			if PP_dist_norm < 2 # particles are overlapping	
				r_c_DEM = sqrt(radius^2 - (PP_dist/2)^2) # Contact radius as calc'd by DEM		
				q_ij_cond = 2*k_p* r_c_DEM * c_factor_PP * (T_curr[j] - T_curr[i]) # PP conduction (Batchelor's equation), with c_factor_PP used to reduce contact area
			else		
				q_ij_cond = 0				
			end



			# PP CONVECTION (aka PFP conduction)			
			if (PP_dist_norm < 3) # only proceed if the h/r < 0.5, or PP_dist < 3 radii, or else ht is negligible
			
				T_eval = (T_curr[i] + T_curr[j]) / 2 # Choose temp to evaluate the PFP heat transfer. Here it is the mean temp of the two particles
				
				# Correct for reduced Young's mod if overlapping

				if PP_dist_norm < 2 # if overlapping
					PP_dist_real = 2*sqrt( radius^2 - c_factor_PP^2 * (radius^2 - (PP_dist/2)^2) )  # calc real distance using the c_factor
				else
					PP_dist_real = PP_dist
				end
				
				# Bilinear interpolation based on distance and temp with function bilin_interp(row values, col values, data, row_i, col_i)
				q_ij_conv = (T_curr[j] - T_curr[i]) * bilin_interp(PFP_dist, PFP_Temps, PFP_htc_values, PP_dist_real, T_eval, "PP convection") 
			else
				q_ij_conv = 0
			end	


			# PP RADIATION			
			# Special case for near-wall PP radiation	
#			if zPos[i] < 5*radius && zPos[j] < 5*radius # if BOTH particles are near the z=0 wall 
#				VlF = VlF_near_wall # Use the near-wall VlF for this PP radiation calc
#			else
#				VlF = VlF_bulk				
#			end

			VlF = VlF_bulk

			# If particles are overlapping, reset to 2
			if PP_dist_norm < 2.0
				PP_dist_norm = 2.0
			end


			# Find the applicable RDF with a bilinear interpolation over PP distance and volume fraction
			PP_RDF = bilin_interp(PP_RDF_dist, PP_RDF_VlF, PP_RDF_values, PP_dist_norm, VlF, "PP radiation")
			
			# Calculate radiative heat transfer
			q_ij_rad =  PP_RDF * e_p_eqn * area * sigma * (T_curr[j]^4 - T_curr_i_4) # [Watts] q is positive if entering particle i, meaning j is hotter than i
			
			
			# Total heat transfer from particle j to i
			q_PP_tot = q_ij_rad + q_ij_cond + q_ij_conv		


			q_PP[i] = q_PP[i] + q_PP_tot # Add this heat transfer to particle i.
			q_PP[j] = q_PP[j] - q_PP_tot # Subtract from particle j
			

		end # only do calcs if the particle pair is closer than PP_cutoff
		
	end # iterate through j particles

end # iterate through i particles

return q_PP # total ht rate to each i particle

end

