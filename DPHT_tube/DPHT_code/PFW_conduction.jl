# Finds the htc at various temperature and PP center-to-center distances. 
# Does numerical integration with the trapezoid rule.  

function PFW_conduction(VlF, k_p, R)

# VlF is volume fraction
# k_p is solid thermal conductivity
# R is particle radius

println("\nExecuting PFW_conduction function.\nExporting htc file to PFP_PFW_conduction folder.")

# Temps and distances to create htc table
Temp_table = [273; 373; 473; 573; 673; 773; 873; 973; 1073; 1173; 1273]
dist_all = collect(0.8*R : 0.02*R : 1.52*R)
n_distances = size(dist_all,1)
htc_values = zeros(  size(dist_all,1), size(Temp_table,1) )

# Interpolate table from Incropera and Dewitt for thermal conductivity of air at Temp
Temp_all = [300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1100, 1200, 1300] # Kelvin
k_air_all = [0.0263, 0.030, 0.0338, 0.0373, 0.0407, 0.0439, 0.0469, 0.0497, 0.0524, 0.0549, 0.0573, 0.0596, 0.062, 0.0643, 0.0667, 0.0715, 0.0763, 0.082] # W/mK
k_f_interp = Spline1D(Temp_all, k_air_all) # Set up 1D interpolation with Dierckx to find K_f at Temp

for temp_count = 1:size(Temp_table,1)

	Temp = Temp_table[temp_count]

	k_f = k_f_interp(Temp)
	
	#println(K_f)


	phi = 1 - VlF # porosity


	q_int_tot = zeros(n_distances,1) # Integrated value, the htc

	# Iterate over particle center-to-center distance
	for i = 1:n_distances
		 
		d_PW = dist_all[i] # PW dist

		r_ij = 0.560*R*(1-phi)^(-1/3)

		H = (d_PW - R) # H is positive if they aren't touching, H is negative if they're overlapping

		H_R_ratio = H/R # When this is above 0.5 it can be neglected (particles are too far apart for meaningful heat transfer)

		# Integral limits, based on d_PW. r_bot is lower integration limit
		if H > 0 # not overlapping
			 r_bot = 0
		else # overlapping			  
			 r_bot = sqrt(R^2 - d_PW^2) #lower integration limit
		end

		r_top = r_ij*R / sqrt( r_ij^2 + (R+H)^2 ) # upper integration limit

		first_iter = true
		delta_r = 0.0000001 # radial increment size, meters
		r = r_bot : delta_r : r_top
		n_evals = size(r,1) # There is one more evaluation than there are trapezoids
	
		q_flux_r = zeros(n_evals,1)
		q_ht_r = zeros(n_evals,1)
		q_int = zeros(n_evals,1)
	
		for n=1:n_evals
	
			# Using the Cheng/Zulli equation, but take off the 2*pi*r*dr to just get the flux as a function of radius, so it can be double checked
#			q_flux_r[n] = 1 / ( ( sqrt(R^2 - r[n]^2) - (r[n]/r_ij)*(R+H))*(1/k_p + 1/k_w) + (2/K_f)*( (R + H) - sqrt(R^2 - r[n]^2))) 
			
			l_solid = sqrt(R^2-r[n]^2) - r[n] * (R+H)/r_ij
			l_fluid = (R+H) - sqrt(R^2 - r[n]^2)
			q_flux_r[n] = 1 / (l_solid/k_p + l_fluid/k_f)
			
			# The function to integrate has 2 pi r dr
			q_ht_r[n] = q_flux_r[n] * 2*pi*r[n]
		
			# Use the trapezoid rule to find the area under the curve. q_int is the area in each trapezoid. 
			# n=1 has nothing in it. n=2 is for the trapezoid between the 1st and 2nd r points, etc. 

		
			if first_iter == false
				 q_int[n] =  ((q_ht_r[n] + q_ht_r[n-1] ) / 2 ) * delta_r     
			end
		
			n=n+1
			first_iter = false
		
		end

		q_int_tot[i] = sum(q_int) # This is the integral from r_bot to r_top, or the htc. It needs to be multiplied by delta T to get actual ht in W.

		htc_values[i,temp_count] = q_int_tot[i]

	end # end iterate over distances

end # end iterate through table temperatures



all_outputs = hcat(dist_all,htc_values)
Temp_table = Temp_table'
temps_row = hcat(0, Temp_table) 
all_outputs = vcat(temps_row,all_outputs)

open("../PFP_PFW_conduction/PFW_htc_values.txt", "w") do io
	for i = 1:size(all_outputs,1)
		for j = 1:size(all_outputs,2)
			@printf(io, "%.8e", all_outputs[i,j])
		
			if j == size(all_outputs,2)
				@printf(io, "\n")
			else
				@printf(io, " ")
			end
		
		end
	end
end




end;
