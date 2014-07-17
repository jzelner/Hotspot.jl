module DBM

using Distributions
using DataFrames

export smoothing_window

# random_hotspot <- function(x,y, h2, in_p, out_p) {
	
# 	h <- h2/2
# 	mid_x <- (max(x) + min(x))/2
# 	mid_y <- (max(y) + min(y))/2

# 	hp_x <- c(mid_x-h, mid_x+h)
# 	hp_y <- c(mid_y-h, mid_y+h)

# 	in_spot <- (x > hp_x[1] & x < hp_x[2]) & (y > hp_y[1] & y < hp_y[2])

# 	z <- rbinom(length(x), 1, out_p)
# 	z[in_spot] <- rbinom(sum(in_spot), 1, in_p)

# 	vals <- list(z = z, coord = c(hp_x, hp_y) )
# 	return(vals)

# }

type Dataset 
	x::Array{Float64}
	y::Array{Float64}
	z::Array{Int}
end

function random_hotspot(x::Array{Float64}, y::Array{Float64}, width::Float64, in_p::Float64, out_p::Float64)
	
	h = width/2.0
	z = Array(Int, length(x))
	in_rand = Bernoulli(in_p)
	out_rand = Bernoulli(out_p)


	mid_x = (maximum(x) + minimum(x)) / 2.0
	mid_y = (maximum(y) + minimum(y)) / 2.0

	hp_x = (mid_x-h, mid_x+h)
	hp_y = (mid_y-h, mid_y+h)

	for i=1:length(x)
		xx = x[i]
		yy = y[i]
		if xx >= hp_x[1] && xx <= hp_x[2] && yy >= hp_y[1] && yy <= hp_y[2]
				z[i] = rand(in_rand)
			else
				z[i] = rand(out_rand)
		end

	end

	return Dataset(x,y,z)
end


function circle_points(x::Array{Float64}, y::Array{Float64}, r::Float64) 

	mid_x = (maximum(x) + minimum(x)) / 2
	mid_y = (maximum(y) + minimum(y)) / 2
	
	offsets = 0:(2*pi/20):(2*pi)-(2*pi/20)
	vals = ["x" => Float64[], "y" => Float64[]]

	for i=offsets
		push!(vals["x"], mid_x + r*sin(i)) 
		push!(vals["y"], mid_x + r*cos(i))
	end

	return vals
end

function all_dist(x::Array{Float64},y::Array{Float64}, cx::Float64, cy::Float64) 
	dist = Array(Float64,length(x))
	for i = 1:length(x)
		dist[i] = sqrt((x[i]-cx)^2 + (y[i]-cy)^2)
	end
	return dist
end



function risk_grid(minX::Float64, maxX::Float64, minY::Float64, maxY::Float64, dim::Int) 
	x = minX:((maxX-minX)/(dim-1)):maxX 
	y = minX:((maxX-minX)/(dim-1)):maxY
	vals = ["x" => Float64[], "y" => Float64[]]

	for i = x
		for j = y
			push!(vals["x"],i)
			push!(vals["y"],j)
		end
	end

	return vals
end

function risk_grid(x::Array{Float64}, y::Array{Float64}, dim::Int)
	risk_grid(minimum(x), maximum(x), minimum(y), maximum(y), dim)
end

function risk_grid(ds::Dataset, dim::Int) 
	risk_grid(ds.x, ds.y, dim)
end

function ecdf(x::Array{Float64}, v::Float64) 
	if !issorted(x)
		sort!(x)
	end
	return searchsortedlast(x,v) / length(x)
end

function ecdf(x::Array{Float64}, v::Array{Float64}) 
	vals = Array(Float64, length(v))
	if !issorted(x)
		sort!(x)
	end

	for i = 1:length(v)
		vals[i] = searchsortedlast(x,v[i]) / length(x)
	end
	return vals
end


function smoothing_window{T<:Number}(ed::Array{T,1}, sd::Array{T,1}, p::Float64)

	#Returns a dictionary with a list of lower bound indices and a list of upper
	#bound indices
	sort!(ed)
	vals = ["lb" => Float64[], "ub" => Float64[], "density" => Float64[]]

	#Calculate the number of cases on either side of the selected case we'll use for 
	#calculating indices

	# print(length(ed))
	stride_length = convert(Int, ceil((p / 2) * length(ed)))
	for d=sd
		i = searchsortedlast(ed, d)
		lb = i-stride_length

		if lb < 1
			lb = 1
		end

		ub = i + stride_length
		if ub > length(ed)
			ub = length(ed)
		end

		# num_ed = length(ed)
		# min_d = minimum(ed)
		# max_d = maximum(ed)
		# print("length = $num_ed, $min_d, $max_d, $stride_length, d = $d, i = $i, LB = $lb, UB = $ub\n")

		push!(vals["lb"], ed[lb])
		push!(vals["ub"], ed[ub])
		push!(vals["density"], (ub-lb)/length(ed))

	end

	return vals

end

function smoothing_window{T<:Number}(ed::Array{T}, sd::Array{T}, p::Float64) 

	num_cp = size(ed)[2]
	num_sd = size(sd)[1]
	num_ed = size(ed)[1]


	lb = Array(Float64, num_sd, num_cp)
	ub = Array(Float64, num_sd, num_cp)
	density = Array(Float64, num_sd, num_cp)
	
	for i=1:num_cp
		sw = smoothing_window(ed[:,i], sd[:,i], p)
		lb[:,i] = sw["lb"]
		ub[:,i] = sw["ub"]
		density[:,i] = sw["density"]
	end

	return ["lb" => lb, "ub" => ub, "density" => density]
end

# function point_score(sw::Dict{ASCIIString, Array{Float64}}, case_ids::Array{Int}) 

# end 

function risk_difference(ddist::Array{Float64}, lb::Array{Float64}, ub::Array{Float64}, density::Array{Float64})

	n_cp = size(ddist)[2]
	n_dp = size(ddist)[1]
	n_sp = size(lb)[1]
	rd = zeros(n_sp)

	for i=1:n_cp

		ddi = ddist[:,i]
		for j=1:n_sp
			ubv = ub[j,i]
			lbv = lb[j,i]
			if ubv > ddi[1]
				if lbv < ddi[1]
					low_index = 1
				else
					low_index = searchsortedlast(ddi, lbv)
				end

				high_index = searchsortedlast(ddi, ubv)
				rd[j] = rd[j] + ((high_index - low_index) / n_dp) - density[j,i]
			end
		end
	end

	return rd
end


function hsmap(ds::Dataset, p::Float64 = 0.1, r::Float64 = 1.0, dim::Int = 100)
	#Get a set of circle points
	cp = circle_points(ds.x, ds.y, r)
	num_cp = length(cp["x"])

	#Make a grid of points with the specified width
	rg = risk_grid(ds, dim)

	num_rg = length(rg["x"])
	rg_x = rg["x"]
	rg_y = rg["y"]

	case_x = ds.x[ds.z .== 1]
	case_y = ds.y[ds.z .== 1] 

	#Calculate the distance of each point on the smoothing grid
	#from each external point

	point_distances = Array(Float64, num_rg, num_cp)
	data_distances = Array(Float64, length(ds.x), num_cp)
	case_distances = Array(Float64, length(case_x), num_cp)
	for i=1:num_cp
		cp_x = cp["x"][i]
		cp_y = cp["y"][i]
		point_distances[:,i] = all_dist(rg_x, rg_y, cp_x, cp_y)
		data_distances[:,i] = sort(all_dist(ds.x, ds.y, cp_x, cp_y))
		case_distances[:,i] = sort(all_dist(case_x, case_y, cp_x, cp_y))
	end

	sw =  smoothing_window(data_distances, point_distances, p)
	rd = risk_difference(case_distances, sw["lb"], sw["ub"], sw["density"])

	min_scores = Float64[]
	max_scores = Float64[]
	for i = 1:100
		#randomize case designations
		case_designations = shuffle(ds.z)
		case_x = ds.x[case_designations .== 1]
		case_y = ds.y[case_designations .== 1] 
		for i=1:num_cp
			case_distances[:,i] = sort(all_dist(case_x, case_y, cp["x"][i], cp["y"][i]))
		end
		
		h = risk_difference(case_distances, sw["lb"], sw["ub"], sw["density"])
		push!(min_scores,minimum(h))
		push!(max_scores,maximum(h))
		
	end

	print(min_scores)

	return DataFrame(x = rg_x, y = rg_y, score = rd), min_scores, max_scores


end


end