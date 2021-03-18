function bilin_interp(x_vals, y_vals, data, xi, yi, calling_function)

# For the matrix "data", the rows are the x_vals and columns are y_vals

#display(x_vals)
#display(y_vals)
#display(data)

x_id_1 = 0 # id of x value in the table that's below xi
x_id_2 = 0 # id of x value in the table that's above xi
y_id_1 = 0
y_id_2 = 0 


for i = 1:(size(x_vals,1)-1) # go from the first to the next-to-last value. 

	if xi >= x_vals[i] && xi <= x_vals[i+1] # Must include the equal to, or else the lowest/highest exact values won't be found. 
		x_id_1 = i
		x_id_2 = i+1
		break
		
	end
	 
end

# Check to see if no value was found. Gives an error if this condition is true. Taking out this warning for speed, but if an index[0] error is thrown, it's probably here.
if x_id_1 == 0
	
	x_min = minimum(x_vals)
	x_max = maximum(x_vals)
	println("In bilin_interp, called by function $calling_function, value to find is outside of table in X (rows) direction. Table has min of $x_min, max of $x_max, and the value to find is $xi ") 
end 	



for j = 1:(size(y_vals,1)-1)
	
	if yi >= y_vals[j] && yi <= y_vals[j+1] 
		y_id_1 = j
		y_id_2 = j+1
		break
	end
	
end

if y_id_1 == 0
	y_min = minimum(y_vals)
	y_max = maximum(y_vals)
	println("In bilin_interp, called by function $calling_function value to find is outside of table in Y (columns) direction. Table has min of $y_min, max of $y_max, and the value to find is $yi ") 	
end


# Find the values at the x and y table points above and below xi and yi
x1 = x_vals[x_id_1]
x2 = x_vals[x_id_2]
y1 = y_vals[y_id_1]
y2 = y_vals[y_id_2]


Q11 = data[x_id_1,y_id_1]
Q12 = data[x_id_1,y_id_2]
Q21 = data[x_id_2,y_id_1]
Q22 = data[x_id_2,y_id_2]


R1 = ((x2 - xi)/(x2 - x1))*Q11 + ((xi - x1)/(x2 - x1))*Q21 # lin interp between x1 and x2 at y1
R2 = ((x2 - xi)/(x2 - x1))*Q12 + ((xi - x1)/(x2 - x1))*Q22 # lin interp between x1 and x2 at y2

P  = ((y2 - yi)/(y2 - y1))*R1  + ((yi - y1)/(y2 - y1))*R2 # lin interp between R1 and R2 at xi

#println("x1: $x1, x2 $x2, y1 $y1, y2 $y2")
#println("R1: $R1, R2 $R2, P $P")
#println("xi, yi, x1, y1, x2, y2: $xi, $yi, $x1, $y1, $x2, $y2") 
#println("Interpolated value: $P")

return P

end
