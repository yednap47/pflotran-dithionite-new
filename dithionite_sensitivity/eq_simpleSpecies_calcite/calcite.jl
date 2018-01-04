import Pflotran
import PyPlot
plt = PyPlot
import DataFrames
df = DataFrames

fname = "calcite-obs-0.tec"
obsname = "obs"

chemnames = Pflotran.getChemNames(fname, obsname)
charges = Pflotran.getObsCharges(fname,obsname)
conc_species = Pflotran.readObsDataset(fname,chemnames)
results_species = df.DataFrame(name=chemnames,conc=conc_species[end,2:end])
df.sort!(results_species, cols = [df.order(Symbol("conc"))],rev=true)

totnames = Pflotran.getChemNamesTotal(fname,obsname)
conc_tots = Pflotran.readObsDataset(fname,totnames)
results_tots = df.DataFrame(name=totnames,conc=conc_tots[end,2:end])

mnrlnames = Pflotran.getChemNamesMnrl(fname,obsname)
conc_mnrl = Pflotran.readObsDataset(fname,mnrlnames)
results_mnrl = df.DataFrame(name=mnrlnames,conc=conc_mnrl[end,2:end])

pH = Pflotran.readObsDataset(fname,["pH"])[2,end]

for totname in totnames
    if contains(totname,"Cr+++") || contains(totname,"Fe+++")
        print("    $(split(split(totname,"[M]")[1],"Total ")[2])")
        @printf "%.4e " 1.0e-20
        print("T")
        println()
    else
        print("    $(split(split(totname,"[M]")[1],"Total ")[2])")
        @printf "%.4e " results_tots[results_tots[Symbol("name")].==totname,Symbol("conc")][1]
        print("T")
        println()
    end
end

println()
for mnrlname in mnrlnames
    if contains(mnrlname,"Cr(OH)3") || contains(mnrlname,"Fe(OH)3")
        print("    $(split(mnrlname,"VF")[1][4:end])")
        @printf "%.4e " 0.0
        print("1.d3")
        println()
    else
        print("    $(split(mnrlname,"VF")[1][4:end])")
        @printf "%.4e " results_mnrl[results_mnrl[Symbol("name")].==mnrlname,Symbol("conc")][1]
        print("1.d3")
        println()
    end
end

# Now compare with minerals_initial
fname2 = "calcite_initial-obs-0.tec"
conc_tots2 = Pflotran.readObsDataset(fname2,totnames)
conc_mnrl2 = Pflotran.readObsDataset(fname2,mnrlnames)

results_mnrl_compare = df.DataFrame(name=mnrlnames,one=conc_mnrl[end,2:end],two=conc_mnrl2[end,2:end], error = conc_mnrl[end,2:end]-conc_mnrl2[end,2:end])
results_tots_compare = df.DataFrame(name=totnames,one=conc_tots[end,2:end], two=conc_tots2[end,2:end], error = conc_tots[end,2:end]-conc_tots2[end,2:end])

pH2 = Pflotran.readObsDataset(fname2,["pH"])[2,end]

# Now check injectant
fname3 = "calcite_injectant-obs-0.tec"
conc_tots3 = Pflotran.readObsDataset(fname3,totnames)
results_tots3 = df.DataFrame(name=totnames,conc=conc_tots3[end,2:end])

results_tots_compare = df.DataFrame(name=totnames,one=conc_tots3[end,2:end])

pH3 = Pflotran.readObsDataset(fname3,["pH"])[2,end]

println()
for totname in totnames
    print("    $(split(split(totname,"[M]")[1],"Total ")[2])")
    @printf "%.4e " results_tots3[results_tots3[Symbol("name")].==totname,Symbol("conc")][1]
    print("T")
    println()
end
