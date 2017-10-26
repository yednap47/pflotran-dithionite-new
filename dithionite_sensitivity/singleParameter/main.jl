import Mads
import ExcelReaders
xlr = ExcelReaders

function runforabit(command, timelimit, pollinterval=1)
    # kills terminal command if time limit is exceeded
    # note, all times are in seconds
    starttime = now()
    process = spawn(command)
    while !process_exited(process) && float(now() - starttime) / 1000 < timelimit
        sleep(pollinterval)
    end
    if !process_exited(process)
        kill(process)
        return false
    else
        return true
    end
end

# User info
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

rundir =  "attempt1"
nstops = 3 # nstops below, nstops above base

# Default parameters for writing MADS file
fname = "../parameters.xlsx"
sheetname = "mads"
jcommand = "read_data.jl"
soltype = "external"
simbasename = "1d-allReactions-10m-uniformVelocity"
startover = "false"

# pflotran information
basedir = "/lclscratch/sach/Programs/pflotran-dithionite-new/dithionite_sensitivity/singleParameter"
pfle = "/lclscratch/sach/Programs/pflotran-dithionite-new/src/pflotran/pflotran"
np = 8
maxruntime = 15 * 60 # seconds from minutes

# write the mads file
f = xlr.openxl(fname)
paraminfo = xlr.readxl(fname, "$(sheetname)!A2:D12")

# make run directory
if !isdir(joinpath(basedir,rundir))
    mkdir(joinpath(basedir,rundir))
end

outfile = open(joinpath(rundir,"$(simbasename).mads"), "w")
println(outfile,"Julia command: $(jcommand)")
println(outfile,"Parameters:")

for i in 1:length(paraminfo[:,1])
    write(outfile, "- log_$(paraminfo[i,1]): ")
    @printf(outfile, "{init: %0.4f, ", paraminfo[i,2])
    @printf(outfile, "min: %0.4f, ", paraminfo[i,3])
    @printf(outfile, "max: %0.4f, ", paraminfo[i,4])
    write(outfile, "type: opt}\n")
    println(outfile, "- $(paraminfo[i,1]): {exp: \"10^log_$(paraminfo[i,1])\"} ")
end

close(outfile)

# read the mads file
madsdata = Mads.loadmadsfile(joinpath(rundir,"$(simbasename).mads"))

# get param names (non-log) and param ranges
logparams_init = Mads.getparamsinit(madsdata)
logparams_min = Mads.getparamsmin(madsdata)
logparams_max = Mads.getparamsmax(madsdata)
logparamkeys = Mads.getparamkeys(madsdata)
paramkeys = map(x->split(x,"log_")[2],logparamkeys)

# sensitivity analysis
for sensparam in sensparams
    # reinitialize paramkeys and params vals
    paramkeys_act = deepcopy(paramkeys)
    logparams_init_act = deepcopy(logparams_init)
    logparams_min_act = deepcopy(logparams_min)
    logparams_max_act = deepcopy(logparams_max)
    
    # make range for sensitivity analysis
    iparamloc = find(x -> x == sensparam,paramkeys)[1]

    # find indices for fast and slow site parameters
    ik_fe2_o2_fast = find(x -> x == "k_fe2_o2_fast",paramkeys)[1]
    ik_fe2_cr6_fast = find(x -> x == "k_fe2_cr6_fast",paramkeys)[1]
    ik_fe2_o2_slow = find(x -> x == "factor_k_fe2_o2_slow",paramkeys)[1]
    ik_fe2_cr6_slow = find(x -> x == "factor_k_fe2_cr6_slow",paramkeys)[1]
    
    # calculate params for slow sites kinetic rate constants
    logparams_init_act[ik_fe2_o2_slow] =  log10(10^logparams_init[ik_fe2_o2_fast]  * 10^logparams_init[ik_fe2_o2_slow])
    logparams_init_act[ik_fe2_cr6_slow] = log10(10^logparams_init[ik_fe2_cr6_fast] * 10^logparams_init[ik_fe2_cr6_slow])
    logparams_min_act[ik_fe2_o2_slow] =   log10(10^logparams_init[ik_fe2_o2_fast]  * 10^logparams_min[ik_fe2_o2_slow])
    logparams_min_act[ik_fe2_cr6_slow] =  log10(10^logparams_init[ik_fe2_cr6_fast] * 10^logparams_min[ik_fe2_cr6_slow])
    logparams_max_act[ik_fe2_o2_slow] =   log10(10^logparams_init[ik_fe2_o2_fast]  * 10^logparams_max[ik_fe2_o2_slow])
    logparams_max_act[ik_fe2_cr6_slow] =  log10(10^logparams_init[ik_fe2_cr6_fast] * 10^logparams_max[ik_fe2_cr6_slow])

    # rename paramkeys for slow sites kinetic rate constants
    paramkeys_act[ik_fe2_o2_slow] = "k_fe2_o2_slow"
    paramkeys_act[ik_fe2_cr6_slow] = "k_fe2_cr6_slow"

    sensvals = Array{Float64}(0)
    sensvals = append!(sensvals,logparams_init_act[iparamloc])

    # append param values below base value
    smallparams = collect(linspace(logparams_min_act[iparamloc],logparams_init_act[iparamloc],nstops+1)[1:end-1])
    largeparams = collect(linspace(logparams_init_act[iparamloc],logparams_max_act[iparamloc],nstops+1)[2:end])
    sensvals = append!(smallparams,sensvals)
    sensvals = append!(sensvals,largeparams)
    sensvals = 10.^sensvals
    @show sensvals 

    # make a matrix [parameters, sensitivity run]
    paramarray = Array{Float64}(length(paramkeys),length(sensvals))
    for i in 1:length(sensvals)
        paramarray[:,i] = 10.^logparams_init_act
        if sensparam == "k_fe2_o2_fast"
            paramarray[iparamloc,i] = sensvals[i]
            paramarray[ik_fe2_o2_slow,i] = 10^logparams_init[ik_fe2_o2_slow] * sensvals[i]
        elseif sensparam == "k_fe2_cr6_fast"
            paramarray[iparamloc,i] = sensvals[i]
            paramarray[ik_fe2_cr6_slow,i] = 10^logparams_init[ik_fe2_cr6_slow] * sensvals[i]
        elseif sensparam == "factor_k_fe2_o2_slow"
            paramarray[iparamloc,i] = sensvals[i]
        elseif sensparam == "factor_k_fe2_cr6_slow"
            paramarray[iparamloc,i] = sensvals[i]
        else
            paramarray[iparamloc,i] = sensvals[i]
        end
    end

    # now lets make the inputfiles
    if !isdir(joinpath(basedir,rundir,sensparam))
        mkdir(joinpath(basedir,rundir,sensparam))
    end

    for i in 1:length(sensvals)
        if !isdir(joinpath(basedir,rundir,sensparam,"run$i"))
            mkdir(joinpath(basedir,rundir,sensparam,"run$i"))
        end
        parameters = Dict(zip(paramkeys_act, paramarray[:,i]))
        parameters["ina"]=2*parameters["is2o4"] # Manually make a parameter, ina = 2*is2o4
        templatefilename = "$simbasename.in.tpl"
        outputfilename = joinpath(basedir,rundir,sensparam,"run$i/$simbasename.in")
        Mads.writeparametersviatemplate(parameters, templatefilename, outputfilename)
        cd(joinpath(basedir,rundir,sensparam,"run$i"))
        try
            println("starting $sensparam run $i")
            runforabit(`mpirun -np $np $pfle -pflotranin $(simbasename).in > $(simbasename).txt`, maxruntime)
        catch
            warn("$(rundir)/$(sensparam)/run$(i) failed")
        end
        cd(basedir)
    end
end
