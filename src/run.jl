include("dbm.jl")

using DBM
# using DataFrames

z = ARGS[length(ARGS)]

df = readtable(z)
print(head(df))
ds = DBM.Dataset(df[:x],df[:y],df[:z])

z,mn,mx = DBM.hsmap(ds, 0.1, 1.0, 250)
minmax = DataFrame(min = mn, max = mx)

writetable("out.csv",z)
writetable("minmax.csv",minmax)

run(`./plot.R`)