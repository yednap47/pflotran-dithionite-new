import(Mads)
using PyPlot
plt = PyPlot
import(JLD)

md = Mads.loadmadsfile("1d-allReactions-10m-uniformVelocity-efast.mads")

info("Forward run check")
tic()
o = Mads.forward(md)
Mads.plotmatches(md, o, r"Cr6_Obs", filename="plot-init-Cr6.png")
toc()
