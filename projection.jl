# ROUTINE FOR PROJECTING ONTO FIRST-STAGE FEASIBLE REGION


# project on to the first-stage feasible region
function projectOntoFeasibleSet(p0_val::Array{Float64},rp_val::Array{Float64},rm_val::Array{Float64},alpha_val::Array{Float64})

	proj_mod = Model(solver=GurobiSolver(gurobi_env,Presolve=0,OutputFlag=0))

	# first-stage variables
	@variable(proj_mod, pmin[i] <= p0[i=1:numGen] <= pmax_FS[i])
	@variable(proj_mod, 0 <= rp[i=1:numGen] <= pmax_FS[i] - pmin[i])
	@variable(proj_mod, 0 <= rm[i=1:numGen] <= pmax_FS[i] - pmin[i])
	@variable(proj_mod, 0 <= alpha[i=1:numGen] <= 1)

	# projection objective function
	@objective(proj_mod, Min, sum( (p0[i] - p0_val[i])^2 + (rp[i] - rp_val[i])^2 + (rm[i] - rm_val[i])^2 + (alpha[i] - alpha_val[i])^2 for i = 1:numGen))

	# first-stage constraints
	@constraint(proj_mod, sum(p0[i] for i=1:numGen) == sum(meanDemand[j] for j=1:numLoads))
	@constraint(proj_mod, [i=1:numRegGen], p0[i] + rp[i] <= pmax[i])
	@constraint(proj_mod, [i=1:numRegGen], p0[i] - rm[i] >= pmin[i])
	@constraint(proj_mod, [i=1:numGen; genGuarRes[i]], alpha[i] >= min_alpha)
	@constraint(proj_mod, [i=1:numRegGen; !genRes[i]], alpha[i] == 0)
	@constraint(proj_mod, sum(alpha[i] for i=1:numGen) == 1)
	
	if(caseID == "no_wind_reserves")
		@constraint(proj_mod, [i=1:numWindGen], alpha[i+numRegGen] == min_alpha/10.0)
	end
	

	solve_status = solve(proj_mod)	
	
	p0_proj = getvalue(p0)
	rp_proj = getvalue(rp)
	rm_proj = getvalue(rm)
	alpha_proj = getvalue(alpha)
	
	# Terminate if the projection model does not solve even to suboptimality 
	if(solve_status != :Optimal && solve_status != :Suboptimal)
		println("p0_val: ",p0_val,"  p0_proj: ",p0_proj)
		println("rp_val: ",rp_val,"  rp_proj: ",rp_proj)
		println("rm_val: ",rm_val,"  rm_proj: ",rm_proj)
		println("alpha_val: ",alpha_val,"  alpha_proj: ",alpha_proj)
		println("solve_status: ",solve_status)
		throw(ErrorException("Projection problem wasn't solved to (sub)optimality!!!"))
	end
	
	return p0_proj, rp_proj, rm_proj, alpha_proj, solve_status
end