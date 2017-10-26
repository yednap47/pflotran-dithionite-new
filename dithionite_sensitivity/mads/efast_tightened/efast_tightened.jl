import(Mads)
using PyPlot
plt = PyPlot
import(JLD)

Mads.quietoff()

Mads.setprocs(10+1)
reload("Mads")

tic()
# Load mads file, set weights
@everywhere md = Mads.loadmadsfile("1d-allReactions-10m-uniformVelocity-tightened-efast.mads")

# Mads.invobsweights!(md, 0.01)
# 
# # Manually set weights so that they are equal (or zero)
# obsnames = Mads.getobskeys(md)
# weights = Mads.getobsweight(md)
# myweight = maximum(weights)
# @show myweight
# for w = find(weights .!= 0.0)
#   md["Observations"][obsnames[w]]["weight"] = myweight
# end

info("Global sensitivity analysis")
efastresult = Mads.efast(md, N=1, seed=2016)
Mads.plotobsSAresults(md, efastresult, 
                      filename="sensitivity_global.png", 
                      xtitle = "x", ytitle = "y")

toc()
# elapsed time: 220263.713274993 seconds
JLD.save("efastresult.jld","dictionary",efastresult)

# ================ PLOT RESULTS WITH PYPLOT ================================== #
coolnames  = JLD.load("../setup/coolnames.jld","dictionary")

obsdict = md["Observations"]
obskeys=Mads.getobskeys(md)
paramkeys = Mads.getparamkeys(md)
nT = length(obskeys)
nP = length(paramkeys)
d = Array{Float64}(2, nT)
i = 1
mes = Array{Float64}(nP, nT)
tes = Array{Float64}(nP, nT)
varx = Array{Float64}(nP, nT)
for obskey in obskeys
    d[1,i] = obsdict[obskey]["time"]
    d[2,i] = haskey(obsdict[obskey], "target") ? obsdict[obskey]["target"] : NaN
    j = 1
    for paramkey in paramkeys
        mes[j,i] = efastresult["mes"][obskey][paramkey]
        tes[j,i] = efastresult["tes"][obskey][paramkey]
        varx[j,i] = efastresult["var"][obskey][paramkey]
        j += 1
    end
    i += 1
end

majorFormatter = matplotlib[:ticker][:FormatStrFormatter]("%0.1e")
mycmap = plt.get_cmap("Paired",15)
linewidth = 1.5
f, ax = plt.subplots(4, 1, sharex=true, figsize=(8,10))

# Spaghetti plots
restartdir = "1d-allReactions-10m-uniformVelocity-tightened-efast_restart"
forwarddict = Dict()
times = Mads.getobstime(md)
fnames = readdir(restartdir)
for i in 1:length(fnames)
    forwarddict["run$(i)"] = JLD.load(joinpath(restartdir,fnames[i]))
end
for i in 1:length(fnames)
    ax[1][:plot](times,forwarddict["run$(i)"]["vecresult"])
end

# Other plots
for i in 1:nP
    ax[2][:plot](d[1,:],tes[i,:],lw=linewidth,c=mycmap(i),label=coolnames[paramkeys[i]])
    ax[3][:plot](d[1,:],mes[i,:],lw=linewidth,c=mycmap(i))
    ax[4][:plot](d[1,:],varx[i,:],lw=linewidth,c=mycmap(i))
end
ax[1][:set_xlim](0,361)
ax[1][:yaxis][:set_major_formatter](majorFormatter)
ax[2][:set_ylim](0,1)
ax[3][:set_ylim](0,1)
ax[4][:yaxis][:set_major_formatter](majorFormatter)
ax[1][:set_ylabel]("Cr(VI) [M]")
ax[2][:set_ylabel]("Total Effect")
ax[3][:set_ylabel]("Main Effect")
ax[4][:set_ylabel]("Output Variance")
ax[4][:set_xlabel]("Time [day]")

plt.tight_layout(pad=1.0,w_pad=1.0,h_pad=1.0)

for i in 1:4
box = ax[i][:get_position]()
ax[i][:set_position]([box[:x0], box[:y0], box[:width] * 0.75, box[:height]])
end

ax[2][:legend](loc=2, bbox_to_anchor=(1.0, 1.0),fontsize=11)

f[:canvas][:draw]() # Update the figure

plt.savefig("sensitivity_global_pyplot.png",dpi=600)
plt.close()
