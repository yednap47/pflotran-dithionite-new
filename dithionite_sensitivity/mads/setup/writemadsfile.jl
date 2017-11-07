import ExcelReaders
xlr = ExcelReaders
import DataFrames

#------------------------------------------------------------------------------
# Initialize
#------------------------------------------------------------------------------
# general info
basedir = "/lclscratch/sach/Programs/pflotran-dithionite-new"
simbasename = "1d-allReactions-10m-uniformVelocity"
templatename =  "1d-allReactions-10m-uniformVelocity"
targetsbasename = "syntheticdata-nn-skip10"

# Default parameters for writing MADS file
paramfilename = "../../parameters.xlsx"
sheetname = "mads"
jcommand = "read_data.jl"
soltype = "external"
startover = "true"

#------------------------------------------------------------------------------
# The code
#------------------------------------------------------------------------------
# write the mads file
f = xlr.openxl(paramfilename)
paraminfo = xlr.readxl(paramfilename, "$(sheetname)!A2:D12")
outfile = open(joinpath(".","$(simbasename).mads"), "w")
println(outfile,"Julia command: $(jcommand)")

println(outfile,"Observations:")
targetsf = open(joinpath(basedir,"dithionite_sensitivity","mads","setup","syntheticdata","$(targetsbasename)-targets.txt"), "r")
targets = readlines(targetsf)
close(targetsf)
for target in targets
     write(outfile, target)
end

println(outfile,"Parameters:")
for i in 1:length(paraminfo[:,1])
    write(outfile, "- log_$(paraminfo[i,1]): ")
    @printf(outfile, "{init: %0.4f, ", paraminfo[i,2])
    @printf(outfile, "min: %0.4f, ", paraminfo[i,3])
    @printf(outfile, "max: %0.4f, ", paraminfo[i,4])
    write(outfile, "type: opt}\n")
        if paraminfo[i,1] != "factor_k_fe2_o2_slow" && paraminfo[i,1] != "factor_k_fe2_cr6_slow"
            println(outfile, "- $(paraminfo[i,1]): {exp: \"10^log_$(paraminfo[i,1])\"}")
        end
end

# Do the tranformations for slow site kinetic constants
write(outfile, "- k_fe2_o2_slow: {exp: \"10^log_k_fe2_o2_fast*10^log_factor_k_fe2_o2_slow\"}\n")
write(outfile, "- k_fe2_cr6_slow: {exp: \"10^log_k_fe2_cr6_fast*10^log_factor_k_fe2_cr6_slow\"}\n")

# Added function for initial [Na+]
# write(outfile, "- ina: {exp: \"2*10^log_is2o4\"}\n")

# Finish writing mads file
println(outfile,"Solution: $(soltype)")
println(outfile,"Templates:")
println(outfile,"- tmp1: {tpl: $(templatename).in.tpl, write: $(templatename).in}")
println(outfile,"Restart: $(startover)")
close(outfile)
println("finished writing mads file")

# # Open the template file and put the new times in
# targettimesf = open(joinpath(basedir,"dithionite_sensitivity","mads","setup","syntheticdata","$(targetsbasename)-times.txt"), "r")
# targettimes = readlines(targettimesf)
# close(targettimesf)
# 
# outfile2 = open(joinpath(".","$(simbasename).in.tpl"))
# tempstrings = readlines(outfile2)
# close(outfile2)
# 
# timeindex = find(map(x->contains(x,"  TIMES d"),tempstrings))
# for ti in timeindex
#     tempstrings[ti] = targettimes[1]
# end
# 
# outfile3 = open(joinpath(".","$(simbasename).in.tpl"), "w")
# for tempstring in tempstrings
#     print(outfile3,tempstring)
# end
# close(outfile3)
