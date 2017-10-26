using HDF5
import sachFun

function readdata_efast(d)
    
    function runforabit(command, waittime=0, pollinterval=1)
        # kills terminal command if time limit is exceeded
        # note, all times are in seconds
        starttime = now()
        cmdproc = spawn(command)
        if waittime > 0
            timedwait(() -> process_exited(cmdproc), waittime)
            if process_running(cmdproc)
                warn("Simulation failed on worker ??")
                kill(cmdproc)
                return false
            end
        else
            wait(cmdproc)
        end
    end
    
    function parseh5!(casetag::AbstractString, results)
        # User Data
        np = 8
        maxruntime = 2.0 * 60.0 * 60.0 # seconds from hours

        pfpath = "/lclscratch/sach/Programs/pflotran-dithionite/src/pflotran"
        masstag = casetag * "-mas.dat"
        otag = casetag * ".out"
        crtag = "Cr6_Obs_t"

        # remove old files
        run(`rm -f $masstag`)
        run(`rm -f $otag`)

        directive = "mpirun -np $(np) " * pfpath * "/pflotran -pflotranin " * casetag * ".in >barf.txt"
        asdf = `bash -c "$directive"`
        runforabit(asdf, maxruntime)

        # # METHOD 1: USE THE H5 FILE
        # myvar = ["Total_CrO4-- [M]"]
        # coord_name = "X"
        # distance = 9 # meters
        # mydata = sachFun.readh5_1D_obs(htag,myvar,coord_name,distance)
        # 
        # # remove t_0 and t_final
        # for i in 2:length(mydata[:,1])-1
        #     results[crtag * "$(i-1)"] = mydata[i,2]
        # end

        # METHOD 2: USE THE MASS BALANCE FILE
        myvar = ["east CrO4-- [mol/d]"]
        mydata = sachFun.readObsDataset(masstag,myvar)
        obstimes = mydata[:,1]
        targets = -mydata[:,2]
        for i in 1:length(obstimes)
            results[crtag * "$(i)"] = targets[i]
        end

        # Check to make sure you have all of the observations
        obskeys = Mads.getobskeys(md)
        if length(results) != length(obskeys)
            error("Number of simulated results does not match number of obskeys!!")
        end
        
        return results
    end

    results = DataStructures.OrderedDict{AbstractString,Any}()
    parseh5!("1d-allReactions-10m-uniformVelocity-efast", results)
    return results
end
