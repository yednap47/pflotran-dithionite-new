import DataStructures
import JLD
using LaTeXStrings

coolnames  = DataStructures.OrderedDict{AbstractString,LaTeXStrings.LaTeXString}(
"log_k_s2o4_disp"=>L"\mathrm{k_{S_2O_4^{-2}-disp}}",
"log_k_s2o4_o2"=>L"\mathrm{k_{S_2O_4^{-2}-O_2(aq)}}",
"log_k_s2o4_fe3"=>L"\mathrm{k_{S_2O_4^{-2}-Fe(OH)_3(s)}}",
"log_sitefrac"=>L"\mathrm{\beta}",
"log_k_fe2_o2_fast"=>L"\mathrm{k_{\equiv Fe(II)-O_2(aq)}}",
"log_factor_k_fe2_o2_slow"=>L"\mathrm{\gamma_{\equiv Fe(II)-O_2(aq)}}",
"log_k_fe2_cr6_fast"=>L"\mathrm{k_{\equiv Fe(II)-CrO_4^{-2}}}",
"log_factor_k_fe2_cr6_slow"=>L"\mathrm{\gamma_{\equiv Fe(II)-CrO_4^{-2}}}",
"log_is2o4"=>L"\mathrm{[S_2O_4^{-2}]}",
"log_ifeoh3"=>L"\mathrm{wt\%\,Fe(OH)_3(s)}",
"log_q" =>L"\mathrm{q}",
"log_aqfrac" =>L"\mathrm{\alpha}",
"log_k_fe2aq_o2" =>L"\mathrm{k_{Fe^{+2}-O_2(aq)}}",
"log_k_fe2aq_cr6" =>L"\mathrm{k_{Fe^{+2}-CrO_4^{-2}}}",
)

JLD.save("coolnames.jld","dictionary",coolnames)
