include("dbm.jl")

using DBM
using DataFrames

df = DBM.random_hotspot(Float64[0.001:0.001:1], Float64[0.001:0.001:1], 0.1, 0.9, 0.001)


case_x = df.x[df.z .== 1]
case_y = df.y[df.z .== 1]

z,mn,mx = @time DBM.hsmap(df, 0.1, 1.0, 100)
writetable("out.csv",z)