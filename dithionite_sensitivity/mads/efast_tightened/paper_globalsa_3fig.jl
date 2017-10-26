import JLD
using PyPlot
plt = PyPlot
import Mads

efastresult  = JLD.load("efastresult.jld","dictionary")
coolnames  = JLD.load("../setup/coolnames.jld","dictionary")
md = Mads.loadmadsfile("1d-allReactions-10m-uniformVelocity-tightened-efast.mads")

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

# ================ PLOT RESULTS WITH PYPLOT ================================== #
majorFormatter = matplotlib[:ticker][:FormatStrFormatter]("%0.1e")
mycmap = plt.get_cmap("Paired",15)
linewidth = 1
mysize = 9

f, ax = plt.subplots(3, 1, sharex=true, figsize=(6.5,5.25))
for i in 1:nP
    ax[1][:plot](d[1,:],tes[i,:],lw=linewidth,c=mycmap(i),label=coolnames[paramkeys[i]])
    ax[2][:plot](d[1,:],mes[i,:],lw=linewidth,c=mycmap(i))
end

# Spaghetti
md = Mads.loadmadsfile("1d-allReactions-10m-uniformVelocity-tightened-efast.mads")
restartdir = "1d-allReactions-10m-uniformVelocity-tightened-efast_restart"

forwarddict = Dict()
times = Mads.getobstime(md)
fnames = readdir(restartdir)

for i in 1:length(fnames)
    forwarddict["run$(i)"] = JLD.load(joinpath(restartdir,fnames[i]))
    ax[3][:plot](times,forwarddict["run$(i)"]["vecresult"],lw=linewidth)
end

ax[1][:set_ylim](0,1)
ax[1][:set_title](L"\mathrm{(A)}", size =mysize)
ax[2][:set_ylim](0,1)
ax[2][:set_title](L"\mathrm{(B)}", size =mysize)
ax[1][:tick_params](labelsize=mysize)
ax[2][:tick_params](labelsize=mysize)
ax[3][:tick_params](labelsize=mysize)
ax[1][:set_ylabel]("Effect",size = mysize)
ax[2][:set_ylabel]("Effect",size = mysize)
ax[3][:set_ylabel]("Cr(VI) efflux [mol/d]",size = mysize)
ax[3][:set_xlabel]("Time [day]",size = mysize)
ax[3][:set_xlim](0,361)
ax[3][:set_title](L"\mathrm{(C)}", size =mysize)
ax[3][:yaxis][:set_major_formatter](majorFormatter)

# ax[1][:yaxis][:set_major_formatter](majorFormatter)

# ax[1][:set_ylabel]("Cr(VI) [M]")

plt.tight_layout(pad=1.0,w_pad=1.0,h_pad=1.0)

for i in 1:3
box = ax[i][:get_position]()
ax[i][:set_position]([box[:x0], box[:y0], box[:width] * 0.8, box[:height]])
end

ax[1][:legend](loc=2, bbox_to_anchor=(1.0, 1.03),fontsize=mysize-2, frameon=false)

f[:canvas][:draw]() # Update the figure

plt.savefig("paper_globalsa_3fig.png",dpi=600)
plt.close()
