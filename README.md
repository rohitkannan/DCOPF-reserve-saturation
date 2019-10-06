Julia codes for the paper "Stochastic DC Optimal Power Flow With Reserve Saturation". Results in the paper were generated using Julia 0.6.2, JuMP 0.18.2, Gurobi 7.5.2 and IPOPT 3.12.8.

6-bus data are in "6bus.jl", and 118-bus data are in "118bus.jl" and the associated files for load and wind correlation data.

Stochastic approximation algorithm for the smooth approximation model is implemented in "stochapprox.jl". SAA for the generator penalty model is implemented in "gp_model.jl", and SAA for the conservative affine policy model is implemented in "cap_model.jl". 
Run the above files to generate results for these models.


Note: Make sure to change the working directory in all of the above files.