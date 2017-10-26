# Id: test.jl, Mon 19 Jun 2017 05:10:27 PM MDT pandeys #
# Created by Sachin Pandey, LANL
# Description: Use information about Cr(III) from h5 file and bound Fe(II) in
#    mass balance file to create spider plot metrics
#------------------------------------------------------------------------------
import Pflotran
import Mads
import PyPlot
plt = PyPlot
using LaTeXStrings
import JLD

function getTotalCr3_fromh5(istop,filename,myvar)
    # calculate moles using m^3/m^3_bulk, dxyz and molar volume
    results = Pflotran.readh5_1D(filename,myvar,coord_name,mytime)
    volumes = dx*dy*dz
    moles = results[:,2].*volumes/MV
    total_moles = sum(moles)
    return total_moles
end

#------------------------------------------------------------------------------
# User info
#------------------------------------------------------------------------------
basedir = "/lclscratch/sach/Programs/pflotran-dithionite-new/dithionite_sensitivity/singleParameter"
rundir = "attempt1"
simbasename = "1d-allReactions-10m-uniformVelocity"
coolnames = JLD.load("../mads/setup/coolnames.jld","dictionary")

sensparams = [
              "k_s2o4_disp",
              "k_s2o4_o2",
              "k_s2o4_fe3",
              "fraction",
              "k_fe2_o2_fast",
              "factor_k_fe2_o2_slow",
              "k_fe2_cr6_fast",
              "factor_k_fe2_cr6_slow",
              "is2o4",
              "ifeoh3",
              "q",
              ]

logsensparams = [
                "log_k_s2o4_disp",
                "log_k_s2o4_o2",
                "log_k_s2o4_fe3",
                "log_fraction",
                "log_k_fe2_o2_fast",
                "log_factor_k_fe2_o2_slow",
                "log_k_fe2_cr6_fast",
                "log_factor_k_fe2_cr6_slow",
                "log_is2o4",
                "log_ifeoh3",
                "log_q",
                ]

nstops = 3 # number of sensitivity runs
mytime = 365
timeunits = "d"
MV = 33.1/(100)^3 # m^3/mol
coord_name = "X"

# Get grid information for calculating volumes
dx = vcat(0.01*ones(200),0.1*ones(40),0.2*ones(20))
# dx = 0.01*ones(100)
dy = 0.1
dz = 0.1

#------------------------------------------------------------------------------
# Get parameter info from madsfile
#------------------------------------------------------------------------------
madsdata = Mads.loadmadsfile(joinpath(rundir,"$(simbasename).mads"))

# get param names (non-log) and param ranges
logparams_init = Mads.getparamsinit(madsdata)
logparams_min = Mads.getparamsmin(madsdata)
logparams_max = Mads.getparamsmax(madsdata)
logparamkeys = Mads.getparamkeys(madsdata)
paramkeys = map(x->split(x,"log_")[2],logparamkeys)

#------------------------------------------------------------------------------
# Build dictionary of results
#------------------------------------------------------------------------------
# loop for parameters we are interested
results = Dict()
for sensparam in sensparams
    @show sensparam
    results[sensparam] = Dict()

    # make range for sensitivity analysis
    iparamloc = find(x -> x == sensparam,paramkeys)[1]

    sensvals = Array{Float64}(0)
    sensvals = append!(sensvals,logparams_init[iparamloc])

    # append param values below base value
    smallparams = collect(linspace(logparams_min[iparamloc],logparams_init[iparamloc],nstops+1)[1:end-1])
    largeparams = collect(linspace(logparams_init[iparamloc],logparams_max[iparamloc],nstops+1)[2:end])
    sensvals = append!(smallparams,sensvals)
    sensvals = append!(sensvals,largeparams)
    sensvals = 10.^sensvals
    results[sensparam]["sensvals"] = sensvals

    # Get total moles of Cr(OH)3 from h5 file
    # Also figure out whether or not the simulation was a completed
    myvar = ["Cr(OH)3(s)_VF"]
    completed = Array{String}(0)
    sensresults = Array{Float64}(0)
    for istop in 1:nstops*2+1
        filename = joinpath(basedir,rundir,sensparam,"run$istop","$(simbasename).h5")
        try
            totalCr3 = getTotalCr3_fromh5(istop,filename,myvar)
            sensresults = append!(sensresults,[totalCr3])
            completed = append!(completed,["yes"])
        catch
            completed = append!(completed,["no"])
        end
        results[sensparam]["success"] = completed
        results[sensparam]["final Cr(III)"] = sensresults
    end

    # Get the MAXIMUM amount of Fe reduced using the mass balance file
    # Get the TOTA amount of s2o4 consumed using the mass balance file
    myvar = ["Global fast_Fe++", "Global slow_Fe++", "east CrO4-- [mol]", "east CrO4-- [mol/d]", "east Cr+++ [mol]"]
    sensresults_fe2 = Array{Float64}(0) # maximum surface bound fe2
    sensresults_cr6 = Array{Float64}(0) # moles of cr6 that reach the outflow
    for istop in 1:nstops*2+1
        if completed[istop] == "yes"
            filename = joinpath(basedir,rundir,sensparam,"run$istop","$(simbasename)-mas.dat")
            obsresults = Pflotran.readObsDataset(filename,myvar,dataframe=true)
            total_sboundfe2 = obsresults[Symbol("Global fast_Fe++")] + obsresults[Symbol("Global slow_Fe++")]
            sensresults_fe2 = append!(sensresults_fe2,[maximum(total_sboundfe2)])
            # find index of desired time in dataframe
            timeindex = find(x -> x == mytime,obsresults[:Time])
            sensresults_cr6 = append!(sensresults_cr6,-(obsresults[Symbol("east CrO4-- [mol]")][timeindex]))
        end
        results[sensparam]["max Fe(II)"] = sensresults_fe2
        results[sensparam]["tot Cr(VI)"] = sensresults_cr6
    end
end

#------------------------------------------------------------------------------
# Spider plots v2
#------------------------------------------------------------------------------
mysize = 9
linewidth = 1
f, ax = plt.subplots(1,2,figsize=(6.5,2.75))
majorFormatter = plt.matplotlib[:ticker][:FormatStrFormatter]("%0.2e")
mycmap = plt.get_cmap("Paired",length(sensparams)+1); extracolor = "0.5"; colorparam = "ifeoh3"
# mycmap = plt.get_cmap("nipy_spectral",length(sensparams)+1); extracolor = "brown"; colorparam = "q"
for sensparam in sensparams
    if sensparam != "d"
    # find the index of sensparam in paramkeys to get base value in params_init
    i = find(x -> x == sensparam,paramkeys)[1]
    basevalue = 10^logparams_init[i]

    # plot loop
    success_summary = results[sensparam]["success"]
    if length(success_summary[success_summary.=="no"])>0
    else
        # # Cr(III)
        # ax[1][:plot](results[sensparam]["sensvals"]/basevalue, results[sensparam]["final Cr(III)"], marker="x", color=mycmap(i), label=sensparam)
        # ax[1][:set_xlim](minimum(results[sensparam]["sensvals"])/basevalue,maximum(results[sensparam]["sensvals"])/basevalue)
        # ax[1][:yaxis][:set_major_formatter](majorFormatter)
        # ax[1][:set_title]("total reduced chromium")
        # ax[1][:set_xscale]("log")
        # ax[1][:set_yscale]("log")
        # ax[1][:set_xlabel]("Fraction of base value")
        # ax[1][:set_ylabel]("total moles Cr(VI) reduced")
        # # ax[1][:legend](loc=0,frameon=false,fontsize=mysize-2)
        # ax[1][:set_xlim](10.0^-1.0,10.0^1.0)
        # # ax[:set_ylim](0.0,4e-2)

        # cumulative Cr(VI) at the outflow
        if sensparam == colorparam
            ax[1][:plot](results[sensparam]["sensvals"]/basevalue, results[sensparam]["tot Cr(VI)"], marker="x", color=extracolor, label=coolnames["log_$(sensparam)"])
        else
            ax[1][:plot](results[sensparam]["sensvals"]/basevalue, results[sensparam]["tot Cr(VI)"], marker="x", color=mycmap(i), label=coolnames["log_$(sensparam)"])
        end
        ax[1][:set_xlim](minimum(results[sensparam]["sensvals"])/basevalue,maximum(results[sensparam]["sensvals"])/basevalue)
        ax[1][:yaxis][:set_major_formatter](majorFormatter)
        ax[1][:set_title](L"\mathrm{(A)}",fontsize = mysize)
        ax[1][:set_xscale]("log")
        ax[1][:set_yscale]("log")
        ax[1][:set_xlabel](L"\mathrm{Fraction\, of\, base\, value}",fontsize = mysize)
        ax[1][:set_ylabel](L"\mathrm{Cr(VI)\, [mol]}",fontsize = mysize)
        ax[1][:set_xlim](10.0^-1.0,10.0^1.0)
        ax[1][:tick_params](labelsize=mysize)
        
        # maximum surface bound Fe(II)
        if sensparam == colorparam
            ax[2][:plot](results[sensparam]["sensvals"]/basevalue, results[sensparam]["max Fe(II)"], marker="x", color=extracolor, label=coolnames["log_$(sensparam)"])
        else
            ax[2][:plot](results[sensparam]["sensvals"]/basevalue, results[sensparam]["max Fe(II)"], marker="x", color=mycmap(i), label=coolnames["log_$(sensparam)"])
        end
        ax[2][:set_xlim](minimum(results[sensparam]["sensvals"])/basevalue,maximum(results[sensparam]["sensvals"])/basevalue)
        ax[2][:yaxis][:set_major_formatter](majorFormatter)
        ax[2][:set_title](L"\mathrm{(B)}",fontsize = mysize)
        ax[2][:set_xscale]("log")
        ax[2][:set_yscale]("log")
        ax[2][:set_xlabel](L"\mathrm{Fraction\, of\, base\, value}",fontsize = mysize)
        ax[2][:set_ylabel](L"\mathrm{\equiv Fe(II)\, [mol]}",fontsize = mysize)
        ax[2][:set_xlim](10.0^-1.0,10.0^1.0)
        ax[2][:set_ylim](10.0^-1.4,10.0^1.0)
        ax[2][:tick_params](labelsize=mysize)
    end
    end
end
box = ax[1][:get_position]()
ax[1][:set_position]([box[:x0]-0.02, box[:y0]+0.09, box[:width] * 0.8, box[:height] * 0.9])

box = ax[2][:get_position]()
ax[2][:set_position]([box[:x0]-0.05, box[:y0]+0.09, box[:width] * 0.8, box[:height] * 0.9])

ax[2][:legend](loc=0,fontsize=mysize-2,loc=2, bbox_to_anchor=(1.1, 1.03), frameon = false)
f[:canvas][:draw]() # Update the figure
plt.savefig("results_full_$(rundir)_$(mytime)$(timeunits).png",dpi=600)
# plt.close()

savedata = Dict()
savedata["results"] = results
savedata["sensparams"] = sensparams
savedata["logsensparams"] = logsensparams
savedata["logparams_init"] = logparams_init
savedata["paramkeys"] = paramkeys
# savedata["coolnames"] = coolnames
# savedata["extracolor"] = extracolor
# savedata["mycmap"] = mycmap

JLD.save("results_sensitivity_singleparam.jld","dictionary",savedata)
