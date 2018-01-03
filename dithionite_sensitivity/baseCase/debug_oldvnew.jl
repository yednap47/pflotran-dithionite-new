import Pflotran
import PyPlot
plt = PyPlot
using LaTeXStrings
import DataFrames

fig = plt.figure("pyplot",figsize=(15,7))

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
mytime = [200.0]
timeunits = "days"
mysize = 11

mylabel = ["$(time) $(timeunits)" for time in mytime]
mycmap = plt.get_cmap("Paired",length(mytime)+1)

majorFormatter = plt.matplotlib[:ticker][:FormatStrFormatter]("%0.2e")
for i in 1:length(myvar)
	plt.subplot(2,4,i)
	for j in 1:length(mytime)
		results = Pflotran.readh5_1D("1d-allReactions-10m-uniformVelocity.h5",[myvar[i]],coord_name,mytime[j])
		plt.plot(results[:,1],(results[:,2]),label=mylabel[j],c=mycmap(j),linewidth=5,c="b")
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

majorFormatter = plt.matplotlib[:ticker][:FormatStrFormatter]("%0.2e")
for i in 1:length(myvar)
	plt.subplot(2,4,i)
	for j in 1:length(mytime)
		results = Pflotran.readh5_1D("1d-allReactions-10m-uniformVelocity-old.h5",[myvar[i]],coord_name,mytime[j])
		plt.plot(results[:,1],(results[:,2]),label=mylabel[j],c=mycmap(j),ls="--",linewidth=5,c="m")
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
plt.savefig("debug1.png",dpi=600)
# plt.close()

plt.figure()
# Plot curve used for mads sensitivity analysis
myvar = ["east CrO4-- [mol/d]"]
redo = Pflotran.readObsDataset("1d-allReactions-10m-uniformVelocity-mas.dat",myvar;dataframe=false)
plt.plot(redo[:,1],-redo[:,2],"b-",label = "new simulation",linewidth=5)
redo2 = Pflotran.readObsDataset("1d-allReactions-10m-uniformVelocity-old-mas.dat",myvar;dataframe=false)
plt.plot(redo2[:,1],-redo2[:,2],"m--",label = "old simulation",linewidth=5)

plt.xlim(0,365)
plt.legend()
plt.xlabel("Time [d]")
plt.ylabel(myvar[1])
plt.tight_layout()
plt.savefig("debug2.png",dpi=600)
# plt.close()

df = DataFrames.DataFrame(tnew=redo[:,1],crnew=-redo[:,2],told=redo2[:,1],crold=-redo2[:,2])
