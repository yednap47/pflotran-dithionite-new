import(Mads)
using PyPlot
plt = PyPlot
import(JLD)

md = Mads.loadmadsfile("1d-allReactions-10m-uniformVelocity-efast.mads")

info("Forward run check")
o = Mads.forward(md)
Mads.plotmatches(md, o, r"Cr6_Obs", filename="plot-init-Cr6.png")

info("Local sensitivity analysis for known parameters")
Mads.localsa(md, filename="sensitivity_local.png", datafiles=true)

# Manually plot local sa results for paper
mysize = 9
coolnames  = JLD.load("../setup/coolnames.jld","dictionary")
eigenmat_raw = readdlm("sensitivity_local-eigenmatrix.dat")

eigenmat_names = eigenmat_raw[:,1]
for i in 1:length(eigenmat_names)
    eigenmat_names[i] = coolnames[eigenmat_names[i]]
end

f, ax = plt.subplots(figsize=(4,3))

eigenmat_clean = eigenmat_raw[:,2:end]
eigenmat_clean = map(x->Float64(x),eigenmat_clean)
plt.pcolor(eigenmat_clean[end:-1:1,1:1:end],cmap="RdBu_r",vmin=-1, vmax=1)
cbar = plt.colorbar()
plt.xlim(0,length(eigenmat_names))
plt.ylim(0,length(eigenmat_names))
plt.xlabel(L"\mathrm{Eigenvectors}",fontsize=mysize)
plt.ylabel(L"\mathrm{Parameters}",fontsize=mysize)
plt.xticks(0.5:1:13.5,(L"1",L"2",L"3",L"4",L"5",L"6",L"7",L"8",L"9",L"10",L"11",L"12",L"13",L"14"))
plt.yticks(0.5:1:13.5,eigenmat_names[end:-1:1])

plt.tick_params(axis="x", which="both", bottom="off", top="off")
plt.tick_params(axis="y", which="both", left="off", right="off")
ax[:tick_params](labelsize=mysize)
cbarvalues = -1.0:.2:1.0
cbarvalues = [LaTeXStrings.LaTeXString("\$$cbarvalue\$") for cbarvalue in cbarvalues]

cbar[:ax][:set_yticklabels](cbarvalues)
cbar[:ax][:tick_params](labelsize=mysize-2)

plt.tight_layout()
plt.draw()
plt.savefig("sensitivity_local-eigenmatrix-sach.png",dpi=600)
plt.close()

Mads.quietoff()

Mads.setprocs(10+1)
reload("Mads")

tic()
# Load mads file, set weights
@everywhere md = Mads.loadmadsfile("1d-allReactions-10m-uniformVelocity-efast.mads")

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
restartdir = "1d-allReactions-10m-uniformVelocity-efast_restart"
forwarddict = Dict()
times = Mads.getobstime(md)
fnames = readdir(restartdir)
fnames = fnames[test.==false]
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
