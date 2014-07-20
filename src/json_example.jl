include("dbm.jl")

using DBM
using JSON
using Distributions

xpoints = rand(5000)
ypoints = rand(5000)

df = DBM.random_hotspot(xpoints, ypoints, 0.1, 0.8, 0.1)

rundata = json(df)
df = DBM.from_json(rundata)

result = DBM.hsmap(df, 0.1, 1.0, 100)

json_result = json(result)

# z = DataFrame(x = result["x"], y = result["y"], score = result["score"])

# writetable("out.csv",z)
# writetable("data.csv", outdf)
# writetable("minmax.csv",minmax)

# run(`./plot.R`)