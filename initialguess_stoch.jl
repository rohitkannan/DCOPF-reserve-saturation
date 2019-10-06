# INITIAL GUESS FOR STOCHASTIC APPROXIMATION ALGORITHM BASED ON A GP MODEL


# returns initial guess given line penalty parameters
function getInitialGuess(gamma_line::Float64,delta_line::Float64,useGurobi::Bool=true)

# set generator bounds violation penalty parameters
	gamma_gen::Float64 = 20.0
	delta_gen::Float64 = 1.0

# generate scenarios
	numScenarios::Int64 = 500
	omega_d, omega_w = generateScenarios(numScenarios)

# compute sum of demand deviations for each scenario
	sumMeanDemands::Float64 = sum(meanDemand[j] for j=1:numLoads)
	sumDemands = (sum(omega_d,1))'
	sumDemandDeviations = sumDemands - sumMeanDemands*ones(Float64,numScenarios)
		
# compute the maximum wind power outputs for the scenarios
	pmax_scen = pmax*ones(Float64,1,numScenarios)
	pmax_scen[numRegGen+1:end,:] = omega_w


	if(useGurobi)
		ig_mod = Model(solver=GurobiSolver(gurobi_env,Presolve=0,OutputFlag=0,BarConvTol=1E-06,OptimalityTol=1E-06))
	else
		ig_mod = Model(solver=IpoptSolver(tol=1E-06,print_level=0))
	end

	
# first-stage variables
	@variable(ig_mod, pmin[i] <= p0_ig[i=1:numGen] <= pmax_FS[i])
	@variable(ig_mod, 0 <= rp_ig[i=1:numGen] <= pmax_FS[i] - pmin[i])
	@variable(ig_mod, 0 <= rm_ig[i=1:numGen] <= pmax_FS[i] - pmin[i])
	@variable(ig_mod, 0 <= alpha_ig[i=1:numGen] <= 1)

# recourse variables
	@variable(ig_mod, p_ig[1:numGen,1:numScenarios])
	@variable(ig_mod, theta_ig[1:numBuses,1:numScenarios])

# auxiliary variables
	@variable(ig_mod, posResViol_ig[1:numGen,1:numScenarios] >= 0)
	@variable(ig_mod, negResViol_ig[1:numGen,1:numScenarios] >= 0)
	@variable(ig_mod, lineViol_ig[1:numLines,1:numScenarios] >= 0)
	@variable(ig_mod, genBoundViol_ig[1:numGen,1:numScenarios] >= 0)


	@objective(ig_mod, Min, sum(genCost[i]*p0_ig[i] + posResCost[i]*rp_ig[i] + negResCost[i]*rm_ig[i] for i=1:numGen) + (1.0/numScenarios)*sum(gamma_res*posResCost[i]*posResViol_ig[i,scen] + gamma_res*negResCost[i]*negResViol_ig[i,scen] for scen=1:numScenarios for i = 1:numGen) + (1.0/numScenarios)*sum(gamma_line*(lineViol_ig[k,scen])^2 for k = 1:numLines for scen = 1:numScenarios) + (1.0/numScenarios)*sum(gamma_gen*(genBoundViol_ig[i,scen])^2 for i = 1:numGen for scen = 1:numScenarios))


# first-stage constraints. same as in the reserve saturation model
	@constraint(ig_mod, sum(p0_ig[i] for i=1:numGen) == sumMeanDemands)
	@constraint(ig_mod, [i=1:numRegGen], p0_ig[i] + rp_ig[i] <= pmax[i])
	@constraint(ig_mod, [i=1:numRegGen], p0_ig[i] - rm_ig[i] >= pmin[i])
	@constraint(ig_mod, [i=1:numGen; genGuarRes[i]], alpha_ig[i] >= min_alpha)
	@constraint(ig_mod, [i=1:numRegGen; !genRes[i]], alpha_ig[i] == 0)
	@constraint(ig_mod, sum(alpha_ig[i] for i=1:numGen) == 1)

	if(caseID == "no_wind_reserves")
		@constraint(ig_mod, [i=1:numWindGen], alpha_ig[i+numRegGen] == min_alpha/10.0)
	end


# constraints defining p: affine recourse policy without worrying about generation limits
	@constraint(ig_mod, [i=1:numGen,scen=1:numScenarios], p_ig[i,scen] == p0_ig[i] + alpha_ig[i]*sumDemandDeviations[scen])

# theta constraints
	@constraint(ig_mod, [i=1:numBuses,scen=1:numScenarios], sum(betaAdjBus[i][k]*(theta_ig[i,scen] - theta_ig[linesAdjBus[i][k],scen]) for k = 1:size(linesAdjBus[i],1)) == sum(p_ig[regGenAtBus[i][k],scen] for k = 1:size(regGenAtBus[i],1)) + sum(p_ig[windGenAtBus[i][k],scen] for k = 1:size(windGenAtBus[i],1)) - sum(omega_d[loadAtBus[i][k],scen] for k = 1:size(loadAtBus[i],1)))
# reference bus
	@constraint(ig_mod, [scen=1:numScenarios], theta_ig[numBuses,scen] == 0)

# line violations
	@constraint(ig_mod, [k=1:numLines,scen=1:numScenarios], lineViol_ig[k,scen] >= beta[k]*(theta_ig[lines[k][1],scen] - theta_ig[lines[k][2],scen]) - delta_line*flow_max[k]) 
	@constraint(ig_mod, [k=1:numLines,scen=1:numScenarios], lineViol_ig[k,scen] >= -beta[k]*(theta_ig[lines[k][1],scen] - theta_ig[lines[k][2],scen]) - delta_line*flow_max[k]) 

# auxiliary constraints for violations of positive and negative reserves
	@constraint(ig_mod, [i=1:numGen,scen=1:numScenarios], posResViol_ig[i,scen] >= p_ig[i,scen] - p0_ig[i] - rp_ig[i])
	@constraint(ig_mod, [i=1:numGen,scen=1:numScenarios], negResViol_ig[i,scen] >= p0_ig[i] - p_ig[i,scen] - rm_ig[i])

# generator bounds violation
	@constraint(ig_mod, [i=1:numGen,scen=1:numScenarios], genBoundViol_ig[i,scen] >= p_ig[i,scen] - (pmin[i] + delta_gen*(pmax_scen[i,scen]-pmin[i])))
	@constraint(ig_mod, [i=1:numGen,scen=1:numScenarios], genBoundViol_ig[i,scen] >= (pmin[i] + (1-delta_gen)*(pmax_scen[i,scen]-pmin[i])) - p_ig[i,scen]) 



	solve_status = solve(ig_mod)	
	
	if(solve_status != :Optimal && solve_status != :Suboptimal)
		if(useGurobi)
			println("GP based IG model not solved to (sub)optimality using Gurobi!!!")
		else
			println("GP based IG model not solved to (sub)optimality using Ipopt!!!")
		end
		println("solve_status: ",solve_status)
	end
	
	
	p0_initguess = getvalue(p0_ig)
	rp_initguess = getvalue(rp_ig)
	rm_initguess = getvalue(rm_ig)
	alpha_initguess = getvalue(alpha_ig)


	return p0_initguess, rp_initguess, rm_initguess, alpha_initguess, solve_status
end