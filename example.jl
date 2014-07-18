include("dbm.jl")

using DBM
using DataFrames
using Distributions

xpoints = rand(5000)
ypoints = rand(5000)

df = DBM.random_hotspot(xpoints, ypoints, 0.1, 0.8, 0.1)
outdf = DataFrame(x = df.x, y = df.y, z = df.z)

case_x = df.x[df.z .== 1]
case_y = df.y[df.z .== 1]

result = DBM.hsmap(df, 0.1, 1.0, 100)
minmax = DataFrame(min = result["min"], max = result["max"])

z = DataFrame(x = result["x"], y = result["y"], score = result["score"])

writetable("out.csv",z)
writetable("data.csv", outdf)
writetable("minmax.csv",minmax)

run(`./plot.R`)