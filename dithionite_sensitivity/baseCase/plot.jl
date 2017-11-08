import Pflotran
import PyPlot
plt = PyPlot
using LaTeXStrings

fig = plt.figure("pyplot",figsize=(15,7))
filename = ARGS[1]
# filename = "1d-allReactions-1m-uniformVelocity"

myvar = [
"pH",
"Total_O2(aq) [M]",
"Total_CrO4-- [M]",
"Fe(OH)3(s)_VF",
"Cr(OH)3(s)_VF",
"slow_Fe++ [mol_m^3]",
"fast_Fe++ [mol_m^3]",
"Total_S2O4-- [M]",
]

coord_name = "X"
mytime = [100.0,200.0,300.0]
timeunits = "days"
mysize = 11

mylabel = ["$(time) $(timeunits)" for time in mytime]
mycmap = plt.get_cmap("Paired",length(mytime)+1)

majorFormatter = plt.matplotlib[:ticker][:FormatStrFormatter]("%0.2e")
for i in 1:length(myvar)
	plt.subplot(2,4,i)
	for j in 1:length(mytime)
		results = Pflotran.readh5_1D("$(filename).h5",[myvar[i]],coord_name,mytime[j])
		plt.plot(results[:,1],(results[:,2]),label=mylabel[j],c=mycmap(j))
		ax = plt.gca()
		ax[:tick_params](labelsize=mysize)
		ax[:set_xlabel](L"\mathrm{Distance\; from\; inlet\;[m]}",size=mysize+2)
		ax[:yaxis][:set_major_formatter](majorFormatter)
		ax[:set_ylabel]("$(myvar[i])",size=mysize+2)
		ax[:set_xlim]([0,maximum(results[:,1])])
		# ax[:set_xlim]([0,10])
		# ax[:set_ylim]([0,maximum(results[:,2])])
	end
end

plt.legend(loc=0,frameon=false,fontsize=mysize)
plt.tight_layout(h_pad=.1)
plt.savefig("$(filename)-test.png",dpi=600)
plt.close()

# Plot curve used for mads sensitivity analysis
masstag = "$(filename)-mas.dat"
myvar = ["east CrO4-- [mol/d]"]
redo = Pflotran.readObsDataset(filename * "-mas.dat",myvar;dataframe=false)
plt.plot(redo[:,1],-redo[:,2],"--",label = "new simulation")
plt.xlim(0,365)
plt.legend()
plt.xlabel("Time [d]")
plt.ylabel(myvar[1])
plt.tight_layout()
plt.savefig("$(filename)-mads-sensitivity.png",dpi=600)
plt.close()

# Calculate parameters used in single parameter sensitivity analysis
# Get the MAXIMUM amount of Fe reduced
# Get the TOTAL amount of Cr(VI) that is untreated
mytime = 365.0
myvar = ["Global fast_Fe++", "Global slow_Fe++", "east CrO4-- [mol]", "east CrO4-- [mol/d]", "east Cr+++ [mol]"]
sensresults_fe2 = Array{Float64}(0) # maximum surface bound fe2
sensresults_cr6 = Array{Float64}(0) # moles of cr6 that reach the outflow
obsresults = Pflotran.readObsDataset(masstag,myvar,dataframe=true)
total_sboundfe2 = obsresults[Symbol("Global fast_Fe++")] + obsresults[Symbol("Global slow_Fe++")]
sensresults_fe2 = append!(sensresults_fe2,[maximum(total_sboundfe2)])

# find index of desired time in dataframe
timeindex = find(x -> x == mytime,obsresults[:Time])
sensresults_cr6 = append!(sensresults_cr6,-(obsresults[Symbol("east CrO4-- [mol]")][timeindex]))
max_Fe2 = sensresults_fe2[1]
tot_Cr6 = sensresults_cr6[1]

print("Total Cr(VI) that leaves domain untreated = $(tot_Cr6)\n")
print("Maximum Fe = $(max_Fe2)\n")
