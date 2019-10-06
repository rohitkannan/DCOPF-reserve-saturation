# GENERATOR PENALTY MODEL FORMULATION THAT PENALIZES VIOLATION OF GENERATOR BOUNDS INSTEAD OF USING RESERVE SATURATION


using Distributions, StatsFuns, StatsBase, JuMP, Gurobi, Ipopt


#*========= MODEL INFORMATION ==========
const model = "118bus"
const modelFile = model * ".jl"
include(modelFile)
#*======================================

include("utilities.jl")
include("modelFns.jl")
include("solveRecourse.jl")
include("evaluateSoln.jl")
include("projection.jl")


#*============ GENERATOR BOUNDS PENALTY ==========

# penalty parameters for generator bounds violation
const gamma_gen_val = logspace(0.0,5.0,16)
if(model == "118bus")
	gamma_gen_val = logspace(0.0,4.0,9)
end
const delta_gen = 1.0


#*============ LINE PENALTY ==========
# set line penalty coefficients to get comparable results
const gamma_line_val = logspace(1.0,5.0,9)


#*============ REPLICATES AND SCENARIOS ==========

# number of replicates for each parameter setting
const numReplicates = 5
const startingRep = 1

if(model == "118bus")
	numReplicates = 1
end

# number of scenarios for the SAA problems
const numScenarios = 500


# set seed for reproducibility
srand(4428)


# generate and store Gurobi environment once and for all
const gurobi_env = Gurobi.Env()




# solve the GP SAA model (convex QP) with penalties on generator bounds violation
function solveGPModel(gamma_gen::Float64,delta_gen::Float64,gamma_line::Float64,delta_line::Float64,
				omega_d::Array{Float64},omega_w::Array{Float64},useGurobi::Bool=true)
	
	# compute sum of demand deviations for each scenario
	sumMeanDemands = sum(meanDemand[j] for j=1:numLoads)
	sumDemands = (sum(omega_d,1))'
	sumDemandDeviations = sumDemands - sumMeanDemands*ones(Float64,numScenarios)
	
	# compute the maximum wind power outputs for the scenarios
	pmax_scen = pmax*ones(Float64,1,numScenarios)
	pmax_scen[numRegGen+1:end,:] = omega_w

	
	# initialize model
	if(useGurobi)
		mod = Model(solver=GurobiSolver(gurobi_env,Presolve=0,OutputFlag=0,BarConvTol=1E-06,OptimalityTol=1E-06))
	else
		mod = Model(solver=IpoptSolver(tol=1E-06,print_level=0))
	end

	
	# first-stage variables
	@variable(mod, pmin[i] <= p0[i=1:numGen] <= pmax_FS[i])
	@variable(mod, 0 <= rp[i=1:numGen] <= pmax_FS[i] - pmin[i])
	@variable(mod, 0 <= rm[i=1:numGen] <= pmax_FS[i] - pmin[i])
	@variable(mod, 0 <= alpha[i=1:numGen] <= 1)

	# recourse variables
	@variable(mod, p[1:numGen,1:numScenarios])
	@variable(mod, theta[1:numBuses,1:numScenarios])

	# auxiliary variables
	@variable(mod, posResViol[1:numGen,1:numScenarios] >= 0)
	@variable(mod, negResViol[1:numGen,1:numScenarios] >= 0)
	@variable(mod, lineViol[1:numLines,1:numScenarios] >= 0)
	@variable(mod, genBoundViol[1:numGen,1:numScenarios] >= 0)

	
	# objective function
	@objective(mod, Min, sum(genCost[i]*p0[i] + posResCost[i]*rp[i] + negResCost[i]*rm[i] for i=1:numGen) + (1.0/numScenarios)*sum(gamma_res*posResCost[i]*posResViol[i,scen] + gamma_res*negResCost[i]*negResViol[i,scen] for scen=1:numScenarios for i = 1:numGen) + (1.0/numScenarios)*sum(gamma_line*(lineViol[k,scen])^2 for k = 1:numLines for scen = 1:numScenarios) + (1.0/numScenarios)*sum(gamma_gen*(genBoundViol[i,scen])^2 for i = 1:numGen for scen = 1:numScenarios))

	
	# first-stage constraints. same as in the reserve saturation model
	@constraint(mod, sum(p0[i] for i=1:numGen) == sumMeanDemands)
	@constraint(mod, [i=1:numRegGen], p0[i] + rp[i] <= pmax[i])
	@constraint(mod, [i=1:numRegGen], p0[i] - rm[i] >= pmin[i])
	@constraint(mod, [i=1:numGen; genGuarRes[i]], alpha[i] >= min_alpha)
	@constraint(mod, [i=1:numRegGen; !genRes[i]], alpha[i] == 0)
	@constraint(mod, sum(alpha[i] for i=1:numGen) == 1)
	
	# for the no wind reserves case, set the wind generator participation factor to be small enough
	if(caseID == "no_wind_reserves")
		@constraint(mod, [i=1:numWindGen], alpha[i+numRegGen] == min_alpha/10.0)
	end

	
	# constraints defining p: affine recourse policy without worrying about generation limits
	@constraint(mod, [i=1:numGen,scen=1:numScenarios], p[i,scen] == p0[i] + alpha[i]*sumDemandDeviations[scen])

	# phase angle constraints
	@constraint(mod, [i=1:numBuses,scen=1:numScenarios], sum(betaAdjBus[i][k]*(theta[i,scen] - theta[linesAdjBus[i][k],scen]) for k = 1:size(linesAdjBus[i],1)) == sum(p[regGenAtBus[i][k],scen] for k = 1:size(regGenAtBus[i],1)) + sum(p[windGenAtBus[i][k],scen] for k = 1:size(windGenAtBus[i],1)) - sum(omega_d[loadAtBus[i][k],scen] for k = 1:size(loadAtBus[i],1)))
	# reference bus
	@constraint(mod, [scen=1:numScenarios], theta[numBuses,scen] == 0)

	# line violations
	@constraint(mod, [k=1:numLines,scen=1:numScenarios], lineViol[k,scen] >= beta[k]*(theta[lines[k][1],scen] - theta[lines[k][2],scen]) - delta_line*flow_max[k]) 
	@constraint(mod, [k=1:numLines,scen=1:numScenarios], lineViol[k,scen] >= -beta[k]*(theta[lines[k][1],scen] - theta[lines[k][2],scen]) - delta_line*flow_max[k]) 

	# auxiliary constraints for violations of positive and negative reserves
	@constraint(mod, [i=1:numGen,scen=1:numScenarios], posResViol[i,scen] >= p[i,scen] - p0[i] - rp[i])
	@constraint(mod, [i=1:numGen,scen=1:numScenarios], negResViol[i,scen] >= p0[i] - p[i,scen] - rm[i])

	# generator bounds violation
	@constraint(mod, [i=1:numGen,scen=1:numScenarios], genBoundViol[i,scen] >= p[i,scen] - (pmin[i] + delta_gen*(pmax_scen[i,scen]-pmin[i])))
	@constraint(mod, [i=1:numGen,scen=1:numScenarios], genBoundViol[i,scen] >= (pmin[i] + (1-delta_gen)*(pmax_scen[i,scen]-pmin[i])) - p[i,scen])


	solve_status = solve(mod)		
	
	if(solve_status != :Optimal && solve_status != :Suboptimal)
		if(useGurobi)
			println("GP model not solved to (sub)optimality using Gurobi!!!")
		else
			println("GP model not solved to (sub)optimality using Ipopt!!!")
		end
		println("solve_status: ",solve_status)
	end

	p0_soln = getvalue(p0)
	rp_soln = getvalue(rp)
	rm_soln = getvalue(rm)
	alpha_soln = getvalue(alpha)	
	
	
	return p0_soln, rp_soln, rm_soln, alpha_soln, solve_status
end



#*=======================================
#*============ RUNS START HERE ==========
#*=======================================

tic()

# number of line flow violation penalty iterations
const numGammaLineIter = size(gamma_line_val,1)
# number of generator bounds violation penalty iterations
const numGammaGenIter = size(gamma_gen_val,1)

# print solutions to the screen at each iteration?
const printIterSoln = false


# directory name for storing results
const baseDirName = "C:/Users/rkannan/Desktop/DCOPF-reserve-saturation/experiments/" * model * "/" * "gp_" * caseID
mkpath(baseDirName)


#*========= PRINT OPTIONS ====================
const storeResults = true

const solnFile = "solution.txt"

const p0File = "p0_soln.txt"
const rpFile = "rp_soln.txt"
const rmFile = "rm_soln.txt"
const alphaFile = "alpha_soln.txt"

const costFile = "expected_cost.txt"
const costPowerFile = "expected_cost_power.txt"
const costReservesFile = "expected_cost_reserves.txt"
const costPosResPenFile = "expected_cost_pos_res_pen.txt"
const costNegResPenFile = "expected_cost_neg_res_pen.txt"

const sumLineViolFile = "sumLineViol.txt"
const lineViolProbFile = "lineViolProb.txt"

const powerGenFile = "expected_power_gen.txt"
const windSpillFile = "expected_wind_spill.txt"

const avgLineViolFile = "avgAbsLineFlowViol.txt"
const genBoundViolProbFile = "genBoundViolProb.txt"
const genViolFile = "genBoundViol.txt"

const timeFile = "solutionTime.txt"


#*========= STORE RESULTS ====================
if(storeResults)

	subDirName = baseDirName * "/"
	mkpath(subDirName)
	
	# write details to text file, including some key details about the test instance
	details_file = subDirName * model * ".txt"
	open(details_file, "w") do f
		write(f,"Model: $model \n")
		write(f,"caseID: $caseID \n")
		write(f,"Number of scenarios (numScenarios): $numScenarios \n\n")
		
		write(f,"======================\n")
		write(f,"Generator bounds violation penalty parameters \n")
		write(f,"----------------------\n")
		write(f,"Penalty factors (gamma_gen): $gamma_gen_val \n")
		write(f,"Beyond what fraction to penalize? (delta_gen): $delta_gen \n\n")
		
		write(f,"======================\n")
		write(f,"Line flow violation penalty parameters \n")
		write(f,"----------------------\n")
		write(f,"Penalty factors (gamma_line): $gamma_line_val \n")
		write(f,"Beyond what fraction to penalize? (delta_line): $delta_line \n\n")
		
		write(f,"======================\n")
		write(f,"Reserves parameters \n")
		write(f,"----------------------\n")
		write(f,"How much costlier are reserves compared to normal generation? (reservesCostFactor): $reservesCostFactor \n")
		write(f,"What is the penalty factor for exceeding reserves? (gamma_res): $gamma_res \n")
		write(f,"How much costlier are wind reserves compared to the cheapest regular generator? (windReservesCostFactor): $windReservesCostFactor \n")
		write(f,"\n\n")
		
		write(f,"**********************\n")
		write(f,"SOME PARAMETERS OF THE INPUT MODEL \n\n")
		write(f,"**********************\n")
		
		write(f,"======================\n")
		write(f,"Generator bounds \n")
		write(f,"----------------------\n")
		write(f,"Regular generator upper limits: $regGenMaxPow \n")
		write(f,"Minimum participation factor: $min_alpha \n\n")
		
		write(f,"======================\n")
		write(f,"Generator cost parameters \n")
		write(f,"----------------------\n")
		write(f,"Regular generator linear cost coefficients: $regGenCost \n")
		write(f,"Regular generator positive reserves cost coefficients: $regGenPosResCost \n")
		write(f,"Regular generator negative reserves cost coefficients: $regGenNegResCost \n")
		write(f,"Wind generator positive reserves cost coefficients: $windGenPosResCost \n")
		write(f,"Wind generator negative reserves cost coefficients: $windGenNegResCost \n\n")
		
		write(f,"======================\n")
		write(f,"Line flow parameters \n")
		write(f,"----------------------\n")
		write(f,"Line flow factor: $lineFlowFactor \n")
		write(f,"Line flow limits: $flow_max \n\n")
		
		write(f,"======================\n")
		write(f,"Demand parameters \n")
		write(f,"----------------------\n")
		write(f,"Demand factor: $demandFactor \n")
		write(f,"Average loads/demands: $meanDemand \n")
		write(f,"Demand stdev factor: $demandStdevFactor \n")
		write(f,"Standard deviation of loads/demands: $demandStdev \n")
		write(f,"Demand covariance matrix: $demandCovMat \n\n")
		
		write(f,"======================\n")
		write(f,"Wind parameters \n")
		write(f,"----------------------\n")
		write(f,"Wind factor: $windFactor \n")
		write(f,"Average wind output: $windPredMean \n")
		write(f,"Wind prediction stdev factor: $windPredStdevFactor \n")
		write(f,"Standard deviation of wind output: $windPredStdev \n")
		write(f,"Do the wind generators provide reserves?: $windGenRes \n")
	end	
	
	# write the line penalty coefficients to a file
	gamma_line_file = subDirName * "gamma_line.txt"
	open(gamma_line_file,"w") do f
		for iter_line = 1:numGammaLineIter
			write(f,"$(gamma_line_val[iter_line])\n")
		end
	end
	
	# write the generator penalty coefficients to a file
	gamma_gen_file = subDirName * "gamma_gen.txt"
	open(gamma_gen_file,"w") do f
		for iter_gen = 1:numGammaGenIter
			write(f,"$(gamma_gen_val[iter_gen])\n")
		end
	end
	
end
#*============================================


# replicate the experiments multiple times, if needed
for rep = 1:numReplicates

	println("\n******************************")
	println("******************************")
	println("REPLICATE ",rep,"/",numReplicates)
	println("******************************")
	
	tic()

	
	#*============ SCENARIO GENERATION ==========
	# for each replicate, generate common set of scenarios that are used across the various penalty settings
	omega_d, omega_w = generateScenarios(numScenarios)

	
# in case we want to solve for more replicates later	
	if(rep >= startingRep)

		# sweep through the line limit penalty parameters
		for iter_line = 1:numGammaLineIter

			#*========= STORE RESULTS ====================
			if(storeResults)

				subDirName2 = baseDirName * "/rep" * string(rep) * "/lp" * string(iter_line) * "/"
				mkpath(subDirName2)
				

				cost_file = subDirName2 * costFile
				open(cost_file, "w") do f
				end

				cost_power_file = subDirName2 * costPowerFile
				open(cost_power_file, "w") do f
				end

				cost_reserves_file = subDirName2 * costReservesFile
				open(cost_reserves_file, "w") do f
				end

				cost_pos_res_file = subDirName2 * costPosResPenFile
				open(cost_pos_res_file, "w") do f
				end

				cost_neg_res_file = subDirName2 * costNegResPenFile
				open(cost_neg_res_file, "w") do f
				end

				sum_viol_file = subDirName2 * sumLineViolFile
				open(sum_viol_file, "w") do f
				end

				viol_prob_file = subDirName2 * lineViolProbFile
				open(viol_prob_file, "w") do f
				end

				power_gen_file = subDirName2 * powerGenFile
				open(power_gen_file, "w") do f
				end

				wind_spill_file = subDirName2 * windSpillFile
				open(wind_spill_file, "w") do f
				end

				avg_viol_file = subDirName2 * avgLineViolFile
				open(avg_viol_file, "w") do f
					write(f,"Average/Maximum line flow violations: \n\n")
				end

				gen_viol_prob_file = subDirName2 * genBoundViolProbFile
				open(gen_viol_prob_file, "w") do f
				end

				gen_viol_file = subDirName2 * genViolFile
				open(gen_viol_file, "w") do f
					write(f,"Probability/Median/Maximum generator bound violations: \n\n")
				end

				soln_file = subDirName2 * solnFile
				open(soln_file, "w") do f
				end

				p0_file = subDirName2 * p0File
				open(p0_file, "w") do f
				end

				rp_file = subDirName2 * rpFile
				open(rp_file, "w") do f
				end

				rm_file = subDirName2 * rmFile
				open(rm_file, "w") do f
				end

				alpha_file = subDirName2 * alphaFile
				open(alpha_file, "w") do f
				end

				time_file = subDirName2 * timeFile
				open(time_file, "w") do f
				end
				
			end
			#*============================================

			println("\n==============================")
			println("==============================")
			println("LINE PENALTY ITER ",iter_line,"/",numGammaLineIter," with penalty ",gamma_line_val[iter_line])
			println("==============================")
			
			tic()
	
			# sweep through the penalty coefficients on generator bounds violation
			for iter_gen = 1:numGammaGenIter

				println("\n*********************")
				println("GEN PENALTY ITER ",iter_gen,"/",numGammaGenIter," with penalty ",gamma_gen_val[iter_gen],"\n")

				tic()

				p0_soln, rp_soln, rm_soln, alpha_soln, solve_status = solveGPModel(gamma_gen_val[iter_gen],delta_gen,gamma_line_val[iter_line],delta_line,omega_d,omega_w)
				
				# if Gurobi did not solve the model as expected, try IPOPT
				if(solve_status != :Optimal && solve_status != :Suboptimal)
					println("--- Gurobi did not converge! Going to try IPOPT ---")
					p0_soln, rp_soln, rm_soln, alpha_soln, solve_status = solveGPModel(gamma_gen_val[iter_gen],delta_gen,gamma_line_val[iter_line],delta_line,omega_d,omega_w,false)
				end
				
				solnTime::Float64 = toq()
				
				if(solve_status != :Optimal && solve_status != :Suboptimal)
					println("*****SKIPPING ITERATION BECAUSE BOTH GUROBI AND IPOPT DID NOT CONVERGE!*****")
				else
					p0_soln, rp_soln, rm_soln, alpha_soln, ~ = projectOntoFeasibleSet(p0_soln, rp_soln, rm_soln, alpha_soln)
			
					if(printIterSoln)
						println("GP solution. Regular generator variables first \n")
						println("gamma_line: ",gamma_line_val[iter_line],"  gamma_gen: ",gamma_gen_val[iter_gen])
						println("p0: ",round.(p0_soln,displayPrecision))
						println("rp: ",round.(rp_soln,displayPrecision))
						println("rm: ",round.(rm_soln,displayPrecision))
						println("alpha: ",round.(alpha_soln,displayPrecision))
						println("")
					end

					tic()
					
					# estimate quality of solution
					cost_soln::Float64, cost_power_soln::Float64, lineViolProb_soln::Float64, sumLineViol_soln::Float64, cost_reserves_soln::Float64, cost_pos_res_soln::Float64, cost_neg_res_soln::Float64, power_gen_soln, wind_spill_soln, ~, avgAbsLineFlowViol_soln, maxAbsLineFlowViol_soln = estimateCostOfSoln(p0_soln, rp_soln, rm_soln,alpha_soln,tau_sat,gamma_res,tau_pos,gamma_line_val[iter_line],delta_line,true,true)
					
					probAnyGenBoundViol_soln::Float64, probGenBoundViol_soln, maxGenBoundViol_soln, medianGenBoundViol_soln = estimateGenBoundViolStats(p0_soln,alpha_soln,true,false)
					
					costTime::Float64 = toq()
			
					totalTimeIter::Float64 = solnTime + costTime
					
					println("ITER TIME: ",round.(totalTimeIter,2),"  SOLN TIME: ",round.(solnTime,2),"  COST TIME: ",round.(costTime,2),"\n")

			
					#*========= STORE RESULTS ==========
					if(storeResults)
						
						open(cost_file, "a") do f
							write(f,"$cost_soln \n")
						end
						
						open(cost_power_file, "a") do f
							write(f,"$cost_power_soln \n")
						end
						
						open(cost_reserves_file, "a") do f
							write(f,"$cost_reserves_soln \n")
						end

						open(cost_pos_res_file, "a") do f
							write(f,"$cost_pos_res_soln \n")
						end

						open(cost_neg_res_file, "a") do f
							write(f,"$cost_neg_res_soln \n")
						end
						
						open(sum_viol_file, "a") do f
							write(f,"$sumLineViol_soln \n")
						end
						
						open(viol_prob_file, "a") do f
							write(f,"$lineViolProb_soln \n")
						end
						
						open(power_gen_file, "a") do f
							for i = 1:numGen
								write(f,"$(power_gen_soln[i]) ")
							end
							write(f,"\n")
						end
						
						open(wind_spill_file, "a") do f
							for i = 1:numWindGen
								write(f,"$(wind_spill_soln[i]) ")
							end
							write(f,"\n")
						end
						
						open(avg_viol_file, "a") do f
							write(f,"ITER NUM. $iter_gen \n")
							for k = 1:numLines
								write(f,"($(lines[k][1]),$(lines[k][2])): $(avgAbsLineFlowViol_soln[k])/$(maxAbsLineFlowViol_soln[k])  \n")
							end
							write(f,"\n")
						end
						
						open(gen_viol_prob_file, "a") do f
							write(f,"$probAnyGenBoundViol_soln \n")
						end
						
						open(gen_viol_file, "a") do f
							write(f,"ITER NUM. $iter_gen \n")
							for i = 1:numGen
								write(f,"Gen #$i: $(probGenBoundViol_soln[i])/$(medianGenBoundViol_soln[i])/$(maxGenBoundViol_soln[i]) \n")
							end
							write(f,"\n")
						end

						open(soln_file, "a") do f
							write(f,"ITER NUM. $iter_gen \n")
							write(f,"Expected cost: $cost_soln \n")
							write(f,"Expected cost without line penalty costs: $cost_power_soln \n")
							write(f,"Expected penalty for exceeding positive reserves: $cost_pos_res_soln \n")
							write(f,"Expected penalty for exceeding negative reserves: $cost_neg_res_soln \n")
							write(f,"Expected sum of line flow limit violations: $sumLineViol_soln \n")
							write(f,"ANY line limit violation prob: $lineViolProb_soln \n")
							write(f,"ANY generator bounds violation prob: $probAnyGenBoundViol_soln \n")
							write(f,"p0: $p0_soln \n")
							write(f,"rp: $rp_soln \n")
							write(f,"rm: $rm_soln \n")
							write(f,"alpha: $alpha_soln \n\n")
							write(f,"Expected power generation: $power_gen_soln \n")
							write(f,"Expected wind spill: $wind_spill_soln \n") 
							write(f,"\n\n")
						end

						open(p0_file, "a") do f
							for i = 1:numGen
								write(f,"$(p0_soln[i]) ")
							end
							write(f,"\n")
						end

						open(rp_file, "a") do f
							for i = 1:numGen
								write(f,"$(rp_soln[i]) ")
							end
							write(f,"\n")
						end

						open(rm_file, "a") do f
							for i = 1:numGen
								write(f,"$(rm_soln[i]) ")
							end
							write(f,"\n")
						end

						open(alpha_file, "a") do f
							for i = 1:numGen
								write(f,"$(alpha_soln[i]) ")
							end
							write(f,"\n")
						end

						open(time_file, "a") do f
							write(f,"$totalTimeIter \n")
						end
					end
					#*==================================
				
				end # solve_status

			end # generator penalty iterations
			
			totalTime_1::Float64 = toq()
			
			println("Total time for line penalty iteration #",iter_line,": ",round.(totalTime_1,2))
		
			#*========= STORE RESULTS ==========
			if(storeResults)
				open(time_file, "a") do f
					write(f,"\n\nTOTAL TIME: $totalTime_1 \n")
				end
			end
			#*==================================
			
		end # line penalty iterations

	end # end replicate num check

	replicateTime::Float64 = toq()
	println("Time for replicate #",rep,":  ",round.(replicateTime,2))

end # replicates


toc() # total time for all runs