import Pflotran
using PyCall
@pyimport numpy.random as nr
import PyPlot
plt = PyPlot

function writemadstargets(basedir::String, myvar::Array{String,1}, timeUnits::String; plotresults::Bool=false, h5::Bool=false)
	
	obsdatafname = joinpath(basedir,"baseCase","$(simbasename)-mas.dat")
	if h5
		# put code for reading h5 file here
	else
		results = Pflotran.readObsDataset(obsdatafname,myvar)
	end
	
	if !isdir(joinpath(".","syntheticdata"))
		mkdir(joinpath(".","syntheticdata"))
	end
	outfile = open(joinpath(".","syntheticdata","$(targetsfname)-nn-skip$(skipfactor)-targets.txt"), "w")

	x = results[:,1]
	y1 = -results[:,2]
	if skipfactor != 0
		obstimes = x[1:skipfactor:end]
		targets = y1[1:skipfactor:end]
	else
		obstimes = x
		targets = y1
	end

	# non-negative values only:
	targets_nn = targets[targets.>0]
	obstimes_nn = obstimes[targets.>0]

	if plotresults
		fig = plt.figure("pyplot_fig",figsize=(15,10))
		p = plt.plot(x,y1,linestyle="-",
				marker=" ",
				label= "simulation",
				color = "blue") # Plot a basic line

		plt.scatter(obstimes_nn,targets_nn, c="b", label = "targets, noisy, nn")

		plt.title(myvar)
		plt.legend()
		plt.savefig(joinpath(".","syntheticdata","$(targetsfname)-nn-skip$(skipfactor).png"))
		plt.close()
	end

	for i in 1:length(obstimes_nn)
		 write(outfile, "- $(caltag)_t$(i): ")
		 @printf(outfile, "{target: %0.5e, ", targets_nn[i])
		 @printf(outfile, "weight: 100.0, ")
		 @printf(outfile, "time: %.2f}\n", obstimes_nn[i])
	end
		
	close(outfile)
	
	# print observation times
	outfile2 = open(joinpath(".","syntheticdata","$(targetsfname)-nn-skip$(skipfactor)-times.txt"), "w")
	write(outfile2, "  TIMES $(timeUnits) ")
	for t in obstimes_nn
		@printf(outfile2,"%.1f ", t)
	end
	println(outfile2)
	# println(outfile2,"ncalibrationTargets = $(length(obstimes_nn))")
	close(outfile2)

end

# general info
basedir = "../.."
simbasename = "1d-allReactions-10m-uniformVelocity"
targetsfname = "syntheticdata"

# User info for making synthetic data
myvar = ["east CrO4-- [mol/d]"]
skipfactor = 0 # resolution of targets vs simulation output
caltag = "Cr6_Obs"
timeUnits = "d"

writemadstargets(basedir,myvar,timeUnits,plotresults=true)
