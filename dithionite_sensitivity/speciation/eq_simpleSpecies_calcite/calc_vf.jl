names = ["Fe(OH)3(s)","Cr(OH)3(s)","Calcite"]

# UNITS = cm^3/mol
mv = Dict()
mv["Fe(OH)3(s)"] = 33.1
mv["Cr(OH)3(s)"] = 34.36
mv["Calcite"] = 36.9340

# UNITS = g/mole
mw = Dict()
mw["Fe(OH)3(s)"] = 103.0181
mw["Cr(OH)3(s)"] = 106.87
mw["Calcite"] = 100.0872

rho_bulk = 1200 # kg/m^3
wt_percent = 5 # %
ssa = 1000

vf = Dict()
for name in names
    if name == "Cr(OH)3(s)" || name == "Fe(OH)3(s)"
        vf[name] = 0.0
        print("    $name ")
        @printf "%.2e " vf[name]
        @printf "%.1e" ssa
        println()
    else
        vf[name] = wt_percent/100*rho_bulk*(mv[name]/100.0^3)/(mw[name]/1000)
        print("    $name ")
        @printf "%.2e " vf[name]
        @printf "%.1e" ssa
        println()
    end
end
