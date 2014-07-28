module Hotspot

using Distributions
using JSON
using ProgressMeter

export smoothing_window
export hsmap

type Dataset
	x::Array{Float64}
	y::Array{Float64}
	z::Array{Int}
end

type MapInput
	p::Float64
	r::Float64
	dim::Int
	ds::Dataset
	xbounds
	ybounds
end

function MapInput(p::Float64, r::Float64, dim::Int, ds::Dataset)
	return(MapInput(p,r,dim,ds, None, None))
end

function json(ds::Dataset)
	return json(["x" => ds.x, "y" => ds.y, "z" => ds.z])
end

function from_json(dsarr::Array{Any,1})
	out_arr = MapInput[]
	for ds=dsarr
		x = convert(Array{Float64}, ds["x"])
		y = convert(Array{Float64}, ds["y"])
		z = convert(Array{Int}, ds["z"])
		dim = convert(Int, ds["dim"])
		p = convert(Float64, ds["p"])
		r = convert(Float64, ds["r"])
		if (ds["xbounds"] != "None") & (ds["ybounds"] != "None")
			xbounds = convert(Array{Float64}, ds["xbounds"])
			ybounds = convert(Array{Float64}, ds["ybounds"])
		else 
			xbounds = None
			ybounds = None
		end

		mi = MapInput(p,r,dim, Dataset(x,y,z), xbounds, ybounds)
		push!(out_arr, mi)
	end

	return(out_arr)
end

function from_json(ds::Dict{String, Any})
	x = convert(Array{Float64}, ds["x"])
	y = convert(Array{Float64}, ds["y"])
	z = convert(Array{Int}, ds["z"])

	dim = convert(Int, ds["dim"])
	p = convert(Float64, ds["p"])
	r = convert(Float64, ds["r"])
	mi = MapInput(p,r,dim, Dataset(x,y,z))
	return mi
end

function from_json(ds::String)
	d = JSON.parse(ds)

	x = convert(Array{Float64}, d["x"])
	y = convert(Array{Float64}, d["y"])
	z = convert(Array{Int}, d["z"])

	return Dataset(x,y,z)
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
		push!(vals["y"], mid_y + r*cos(i))
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
	#bound indices and the density of cases in each window
	sort!(ed)
	vals = ["lb" => Float64[], "ub" => Float64[], "density" => Float64[]]

	#Calculate the number of cases on either side of the selected case we'll use for
	#calculating indices

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

		push!(vals["lb"], ed[lb])
		push!(vals["ub"], ed[ub])
		push!(vals["density"], (ub-lb+1)/length(ed))

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
			#If the upper and lower distance bounds for the
			#point are inside the range in where there are actually
			#cases, then we estimate the proportion of cases inside of this window.
			#Otherwise, it's zero and we just go ahead and subtract the control density. 
			if (ubv > ddi[1]) & (lbv < ddi[length(ddi)])
				if lbv <= ddi[1]
					low_index = 1
				else
					low_index = searchsortedlast(ddi, lbv)
				end

				high_index = searchsortedlast(ddi, ubv)

				rd[j] += ((high_index - low_index+1) / length(ddi)) - density[j,i]
			else
				rd[j] -= density[j,i]
			end
		end

	end
	return rd ./ n_cp
end

function hsmap_from_json(fname::String)
	df = Hotspot.from_json(JSON.parsefile(fname))
	result = hsmap(df)
	return JSON.json(result)
end

function hsmap(mi::MapInput) 
	print("Making map with p = $(mi.p), r = $(mi.r), dim = $(mi.dim)\n")
	return hsmap(mi.ds, mi.p, mi.r, mi.dim, mi.xbounds, mi.ybounds)
end

function hsmap(dslist::Array{MapInput})
	returnmaps = Dict{ASCIIString,Array{Float64}}[]
	for ds=dslist
		push!(returnmaps, hsmap(ds))
	end
	return returnmaps
end


function hsmap(ds::Dataset, p::Float64 = 0.1, r::Float64 = 1.0, dim::Int = 100, xbounds = None, ybounds = None)
	#Get a set of circle points
	cp = circle_points(ds.x, ds.y, r)
	num_cp = length(cp["x"])

	#Make a grid of points with the specified width
	if (xbounds != None) & (ybounds != None)
		rg = risk_grid(xbounds[1],xbounds[2],ybounds[1],ybounds[2],dim)
	else
		rg = risk_grid(ds, dim)
	end

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

	nsamp = 100
	min_scores = Array(Float64, nsamp)
	max_scores = Array(Float64, nsamp)
	p = Progress(nsamp, 1, "Taking $nsamp color samples...",convert(Int,ceil(nsamp/10)))
	for i = 1:nsamp
		#randomize case designations
		case_designations = shuffle(ds.z)
		case_x = ds.x[case_designations .== 1]
		case_y = ds.y[case_designations .== 1]
		for j=1:num_cp
			case_distances[:,j] = sort(all_dist(case_x, case_y, cp["x"][j], cp["y"][j]))
		end

		h = risk_difference(case_distances, sw["lb"], sw["ub"], sw["density"])

		min_scores[i] = minimum(h)
		max_scores[i] = maximum(h)
		next!(p)

	end

	return ["x" => rg_x, "y" => rg_y, "score" => rd, "min" => min_scores, "max" => max_scores]

end


end
