# STOCHASTIC APPROXIMATION ALGORITHM FOR SOLVING DCOPF WITH RESERVE SATURATION

using Distributions, StatsFuns, StatsBase, JuMP, Gurobi


# generate and store Gurobi environment once and for all
const gurobi_env = Gurobi.Env()

tic() # overall time

#*========= MODEL INFORMATION ==========
const model = "6bus"
const modelFile = model * ".jl"
include(modelFile)
#*======================================

#*========= INCLUDE SUPPORTING FILES ==============
include("utilities.jl")
include("modelFns.jl")
include("solveRecourse.jl")
include("evaluateSoln.jl")
include("projection.jl")
include("stochgrad.jl")
include("stepLength.jl")
include("initialguess_stoch.jl")


#*============ LINE PENALTY ==========
# set line penalty coefficients to get comparable results
const gamma_line_val = logspace(1.0,5.0,17)
# finer discretization for the no wind reserves case to try and get better approximation of EF
if(caseID == "no_wind_reserves")
	gamma_line_val1 = logspace(1.0,3.0,9)
	gamma_line_val2 = logspace(3.25,3.5,6)
	gamma_line_val3 = logspace(3.75,5.0,6)
	gamma_line_val = [gamma_line_val1;gamma_line_val2;gamma_line_val3]
end

# set seed for reproducibility
srand(1135)


# directory name for storing results
const baseDirName = "C:/Users/rkannan/Desktop/DCOPF-reserve-saturation/experiments/" * model * "/" * "stochapprox_" * caseID
mkpath(baseDirName)


#*========= PRINT OPTIONS ====================
const storeResults = true
const printIterationStats = false
const printFinalSolution = false

const initFile = "initial_guess.txt"
const solnFile = "solution.txt"

const costFile = "expected_cost.txt"
const costPowerFile = "expected_cost_power.txt"
const costReservesFile = "expected_cost_reserves.txt"
const sumLineViolFile = "sumLineViol.txt"
const lineViolProbFile = "lineViolProb.txt"

const avgLineViolFile = "avgAbsLineFlowViol.txt"
const timeFile = "solutionTime.txt"

const bestSolnFile = "best_solution.txt"
const p0BestSolnFile = "p0_best_solution.txt"
const rpBestSolnFile = "rp_best_solution.txt"
const rmBestSolnFile = "rm_best_solution.txt"
const alphaBestSolnFile = "alpha_best_solution.txt"
const bestCostFile = "best_expected_cost.txt"
const bestCostPowerFile = "best_expected_cost_power.txt"
const bestSumLineViolFile = "best_sumLineViol.txt"
const bestLineViolProbFile = "best_lineViolProb.txt"

const bestCostReservesFile = "best_expected_cost_reserves.txt"
const bestCostPosResPenFile = "best_expected_cost_pos_res_pen.txt"
const bestCostNegResPenFile = "best_expected_cost_neg_res_pen.txt"
const powerGenSolnFile = "best_expected_power_gen.txt"
const windSpillSolnFile = "best_expected_wind_spill.txt"


#*========== ALGORITHMIC PARAMETERS ==========
# mini-batch size (m)
const minibatchSize = 20

# max. number of iterations for each replicate
const maxNumIterations = 1000 

# minimum and maximum number of replicates (S)
const minNumReplicates = 20
const maxNumReplicates = 1000 

# number of trials of the algorithm
const numTrials = 5
if(model == "118bus")
	numTrials = 1
end

# improvement requirement for considering another replicate
const impReqForTerm = 1E-05
const termCheckPeriod = 10

# step size modification parameters
const stepLengthDecrCheckIter = 3
const stepLengthIncrCheckIter = 3
const impReqForIncrStepSize = 1E-03
const impReqForDecrStepSize = 1E-03
const stepLengthIncrFactor = 2.0
const stepLengthDecrFactor = 2.0

#*===========================================

# number of line penalty parameter iterations
const numGammaLineIter = size(gamma_line_val,1)


println("")
println("***** STARTING THE STOCHASTIC APPROXIMATION METHOD FOR MODEL ",model," AND CASE ",caseID)
println("")


# run trials/replications
for trial = 1:numTrials

	println("\n******************************")
	println("******************************")
	println("TRIAL ",trial,"/",numTrials)
	println("******************************")

	tic() # trial time


	# sweep through the line limit penalty parameters
	for iter_line = 1:numGammaLineIter

		tic() # line penalty iteration time
		
		gamma_line_iter::Float64 = gamma_line_val[iter_line]

		println("\n==============================")
		println("==============================")
		println("LINE PENALTY ITER ",iter_line,"/",numGammaLineIter," with penalty ",gamma_line_iter)
		println("==============================")

		tic() # preprocessing time
		
		tic()

		#*========== COMPUTE INITIAL GUESS ===========
		cost_init2::Float64 = Inf
		cost_power_init2::Float64 = Inf
		probAnyLineFlowViol_init2::Float64 = 1.0
		sumLineFlowViol_init2::Float64 = Inf
		cost_reserves_init2::Float64 = Inf
		
		# in subsequent iterations, try and use the best solution from the previous iteration
		if(iter_line > 1)
			# project initial guess onto feasible set
			p0_initguess3, rp_initguess3, rm_initguess3, alpha_initguess3, ~ = projectOntoFeasibleSet(p0_init,rp_init,rm_init,alpha_init)

			cost_init2, cost_power_init2, probAnyLineFlowViol_init2, sumLineFlowViol_init2, cost_reserves_init2, ~, ~, ~, ~, ~, ~, ~ = estimateCostOfSoln(p0_initguess3,rp_initguess3,rm_initguess3,alpha_initguess3,tau_sat,gamma_res,tau_pos,gamma_line_iter,delta_line,false,false)
		end

		# compute initial guess for the first-stage variables using the GP model
		p0_initguess, rp_initguess, rm_initguess, alpha_initguess, ig_solve_status = getInitialGuess(gamma_line_iter,delta_line)
				
		if(ig_solve_status != :Optimal && ig_solve_status != :Suboptimal)
			println("--- Gurobi did not solve IG model to optimality! Going to try IPOPT ---")
			p0_initguess, rp_initguess, rm_initguess, alpha_initguess, ig_solve_status = getInitialGuess(gamma_line_iter,delta_line,false)
			if(ig_solve_status != :Optimal && ig_solve_status != :Suboptimal)
				println("--- Ipopt also did not solve IG model to optimality! Going to use a trivial initial guess ---")
				p0_initguess = zeros(Float64,numGen)
				rp_initguess = zeros(Float64,numGen)
				rm_initguess = zeros(Float64,numGen)
				alpha_initguess = ones(Float64,numGen)/numGen
			end
		end

		# project initial guess onto feasible set
		p0_initguess2, rp_initguess2, rm_initguess2, alpha_initguess2, ~ = projectOntoFeasibleSet(p0_initguess,rp_initguess,rm_initguess,alpha_initguess)


		#*========== EVALUATE QUALITY OF INITIAL GUESS ===========
		
		IGSolnTime::Float64 = toq()
		
		tic() # cost time

		cost_init::Float64, cost_power_init::Float64, probAnyLineFlowViol_init::Float64, sumLineFlowViol_init::Float64, cost_reserves_init::Float64, ~, ~, ~, ~, ~, ~, ~ = estimateCostOfSoln(p0_initguess2,rp_initguess2,rm_initguess2,alpha_initguess2,tau_sat,gamma_res,tau_pos,gamma_line_iter,delta_line,false,false)
		
		if(cost_init < cost_init2)
			for i = 1:numGen
				p0_init[i] = p0_initguess2[i]
				rp_init[i] = rp_initguess2[i]
				rm_init[i] = rm_initguess2[i]
				alpha_init[i] = alpha_initguess2[i]
			end
		else
			for i = 1:numGen
				p0_init[i] = p0_initguess3[i]
				rp_init[i] = rp_initguess3[i]
				rm_init[i] = rm_initguess3[i]
				alpha_init[i] = alpha_initguess3[i]
			end
			cost_init = cost_init2
			cost_power_init = cost_power_init2
			probAnyLineFlowViol_init = probAnyLineFlowViol_init2
			sumLineFlowViol_init = sumLineFlowViol_init2
			cost_reserves_init = cost_reserves_init2
		end
		
		println("\n*** INITIAL GUESS ***")
		
		println("Expected cost (incl. all penalties): ",round.(cost_init,displayPrecision))
		println("Expected cost without line penalty costs: ",round.(cost_power_init,displayPrecision))
		println("Expected cost of reserves (incl. reserve penalties): ",round.(cost_reserves_init,displayPrecision))
		println("Expected sum of line flow violations: ",round.(sumLineFlowViol_init,displayPrecision))
		println("Probability that ANY line flow is violated: ",round.(probAnyLineFlowViol_init,displayPrecision),"\n")
		

		costTime::Float64 = toq()
		println("Time for evaluating initial guess: ",round.(IGSolnTime,2),"  time for evaluating quality: ",round.(costTime,2))


		tic() # step length time

		#*============ ESTIMATE STEP LENGTH ==========

		numSamplesForLipEst::Int64 = 200

		stepLength = estimateStepLength(p0_initguess2,rp_initguess2,rm_initguess2,alpha_initguess2,minibatchSize,numSamplesForLipEst,tau_sat,gamma_res,tau_pos,gamma_line_iter,delta_line)

		initialStepLength = 1.0*stepLength

		#*============================================

		stepLengthTime = toq()
		preprocessingTime = toq()

		println("Step length: ",stepLength,"  time for estimating it: ",round.(stepLengthTime,2))
		println("Total preprocessing time: ",round.(preprocessingTime,2),"\n")


		#*========= STORE RESULTS ====================
		if(storeResults)

			subDirName = baseDirName * "/" * "/rep" * string(trial) * "/lp" * string(iter_line) * "/"
			mkpath(subDirName)
			
			# write details to text file, including some key details about the test instance
			details_file = subDirName * model * ".txt"
			open(details_file, "w") do f
				write(f,"Model: $model \n")
				write(f,"caseID: $caseID \n")
				write(f,"Mini-batch size (m): $minibatchSize \n")
				write(f,"Max. number of iterations (N): $maxNumIterations \n")
				write(f,"Min./Max. number of replicates: $minNumReplicates/$maxNumReplicates \n")
				write(f,"improvement requirement for termination: $impReqForTerm \n")
				write(f,"improvement requirement for increasing stepLength: $impReqForIncrStepSize \n")
				write(f,"non-improvement requirement for decreasing stepLength: $impReqForDecrStepSize \n")
				write(f,"check for decreasing step length every $stepLengthDecrCheckIter iterations \n")
				write(f,"check for increasing step length every $stepLengthIncrCheckIter iterations \n")
				write(f,"initial step length: $stepLength \n")
				write(f,"step length increase factor: $stepLengthIncrFactor \n")
				write(f,"step length decrease factor: $stepLengthDecrFactor \n")
				write(f,"number of samples for estimating Lipschitz constant and variance: $numSamplesForLipEst \n")
				write(f,"Preprocessing time: $preprocessingTime \n\n")
				
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
				write(f,"Regular generator upper limits: $regGenMaxPow \n\n")
				
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
				write(f,"Average loads/demands: $meanDemand \n")
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

			# write data to text file
			init_file = subDirName * initFile
			open(init_file, "w") do f
				write(f,"NOTE: regular generator variables come first \n\n")
				write(f,"Expected cost: $cost_init \n")
				write(f,"Expected cost of generation (without line penalty costs): $cost_power_init \n")
				write(f,"Expected sum of line flow limit violations: $sumLineFlowViol_init \n")
				write(f,"Probability that ANY line limit is violated: $probAnyLineFlowViol_init \n\n")
				write(f,"p0: $p0_init \n")
				write(f,"rp: $rp_init \n")
				write(f,"rm: $rm_init \n")
				write(f,"alpha: $alpha_init \n\n")
				write(f,"Time for computing initial guess: ",IGSolnTime)
				write(f,"Time for evaluating quality of initial guess: ",costTime)
			end
	
			# write the line penalty coefficients to a file
			gamma_line_file = subDirName * "gamma_line.txt"
			open(gamma_line_file,"w") do f
				for iter_line2 = 1:numGammaLineIter
					write(f,"$(gamma_line_val[iter_line2])\n")
				end
			end

			# write data to text file
			cost_file = subDirName * costFile
			open(cost_file, "w") do f
			end

			cost_power_file = subDirName * costPowerFile
			open(cost_power_file, "w") do f
			end

			cost_reserves_file = subDirName * costReservesFile
			open(cost_reserves_file, "w") do f
			end

			sum_viol_file = subDirName * sumLineViolFile
			open(sum_viol_file, "w") do f
			end

			viol_prob_file = subDirName * lineViolProbFile
			open(viol_prob_file, "w") do f
			end

			best_cost_file = subDirName * bestCostFile
			open(best_cost_file, "w") do f
			end

			# write data to text file
			best_cost_power_file = subDirName * bestCostPowerFile
			open(best_cost_power_file, "w") do f
			end

			best_sum_viol_file = subDirName * bestSumLineViolFile
			open(best_sum_viol_file, "w") do f
			end

			best_viol_prob_file = subDirName * bestLineViolProbFile
			open(best_viol_prob_file, "w") do f
			end

			avg_viol_file = subDirName * avgLineViolFile
			open(avg_viol_file, "w") do f
			end

			soln_file = subDirName * solnFile
			open(soln_file, "w") do f
			end

			best_soln_file = subDirName * bestSolnFile
			open(best_soln_file, "w") do f
			end

			p0_best_soln_file = subDirName * p0BestSolnFile
			open(p0_best_soln_file, "w") do f
			end

			rp_best_soln_file = subDirName * rpBestSolnFile
			open(rp_best_soln_file, "w") do f
			end

			rm_best_soln_file = subDirName * rmBestSolnFile
			open(rm_best_soln_file, "w") do f
			end

			alpha_best_soln_file = subDirName * alphaBestSolnFile
			open(alpha_best_soln_file, "w") do f
			end

			best_cost_reserves_file = subDirName * bestCostReservesFile
			open(best_cost_reserves_file, "w") do f
			end

			best_cost_pos_res_file = subDirName * bestCostPosResPenFile
			open(best_cost_pos_res_file, "w") do f
			end

			best_cost_neg_res_file = subDirName * bestCostNegResPenFile
			open(best_cost_neg_res_file, "w") do f
			end

			power_gen_soln_file = subDirName * powerGenSolnFile
			open(power_gen_soln_file, "w") do f
			end

			wind_spill_soln_file = subDirName * windSpillSolnFile
			open(wind_spill_soln_file, "w") do f
			end

			time_file = subDirName * timeFile
			open(time_file, "w") do f
			end
			
		end
		#*============================================

		#*========= STORE BEST FOUND SOLUTION ========

		# store the best found solution in terms of expected cost
		# initialize using the initial guess solution
		best_iter::Int64 = 0
		best_cost::Float64 = cost_init
		best_cost_power::Float64 = cost_power_init
		best_lineViolProb::Float64 = probAnyLineFlowViol_init
		best_sumLineViol::Float64 = sumLineFlowViol_init

		p0_best_solution = zeros(Float64,numGen) 
		rp_best_solution = zeros(Float64,numGen) 
		rm_best_solution = zeros(Float64,numGen) 
		alpha_best_solution = zeros(Float64,numGen)
		 
		for i = 1:numGen
			p0_best_solution[i] = p0_init[i]
			rp_best_solution[i] = rp_init[i]
			rm_best_solution[i] = rm_init[i]
			alpha_best_solution[i] = alpha_init[i]
		end

		# store the final iterate and objective for each replicate
		cost_soln = Inf*ones(Float64,maxNumReplicates)
		cost_power_soln = Inf*ones(Float64,maxNumReplicates)
		anyLineViolProb_soln = Inf*ones(Float64,maxNumReplicates)
		sumLineViol_soln = Inf*ones(Float64,maxNumReplicates)

		p0_solution = zeros(Float64,maxNumReplicates,numGen)
		rp_solution = zeros(Float64,maxNumReplicates,numGen)
		rm_solution = zeros(Float64,maxNumReplicates,numGen)
		alpha_solution = zeros(Float64,maxNumReplicates,numGen) 

			
		numStepSizeIncr::Int64 = 0
		numStepSizeDecr::Int64 = 0

		local_best_cost::Float64 = best_cost
		local_best_cost2::Float64 = best_cost
		local_best_cost3::Float64 = best_cost
		
		actualNumIter::Int64 = 0
	

		#*========= OPTIMIZATION LOOP ==========
		for iter_S = 1:maxNumReplicates

			tic()
			objTime::Float64 = 0
			gradTime::Float64 = 0
			projTime::Float64 = 0

			if(printIterationStats)
				println("=====================================")
				println("Optimization iteration # ", iter_S)
			end
			
			# determine actual number of SGD iterations
			uniformRand = rand(1)
			numIterations::Int64 = ceil.(uniformRand[1]*maxNumIterations)
			
			if(printIterationStats)
				println("Number of iterations: ",numIterations)
			end
	
		
			#*========= SET INITIAL GUESS ==========
			p0 = zeros(Float64,numGen)
			rp = zeros(Float64,numGen)
			rm = zeros(Float64,numGen)
			alpha = zeros(Float64,numGen)
			
			# if first replicate, use user-specified initial guess; otherwise, use solution from previous iteration
			if(iter_S == 1)
				for i = 1:numGen
					p0[i] = p0_init[i]
					rp[i] = rp_init[i]
					rm[i] = rm_init[i]
					alpha[i] = alpha_init[i]
				end
			else
				for i = 1:numGen
					p0[i] = p0_solution[iter_S-1,i]
					rp[i] = rp_solution[iter_S-1,i]
					rm[i] = rm_solution[iter_S-1,i]
					alpha[i] = alpha_solution[iter_S-1,i]
				end
			end
			#*======================================

			#*===========================================================
			#* CORE STOCHASTIC APPROXIMATION ALGORITHM IN THIS BLOCK
			#*===========================================================
			for t = 1:numIterations
			
				tic()
			
				#*==== COMPUTE STOCHASTIC GRADIENT =====
				grad_p0, grad_rp, grad_rm, grad_alpha = getStochasticGradient(p0,rp,rm,alpha,minibatchSize,tau_sat,gamma_res,tau_pos,gamma_line_iter,delta_line)
				#*======================================
				
				gradTime += toq()
				
				#*==== SGD STEP =====
				p0 = p0 - stepLength[1]*grad_p0
				rp = rp - stepLength[2]*grad_rp
				rm = rm - stepLength[3]*grad_rm
				alpha = alpha - stepLength[4]*grad_alpha
				#*===================
			
				tic()
				
				#*==== PROJECTION STEP =====
				p0_2, rp_2, rm_2, alpha_2, ~ = projectOntoFeasibleSet(p0,rp,rm,alpha)
				#*==========================
				
				for i = 1:numGen
					p0[i] = p0_2[i]
					rp[i] = rp_2[i]
					rm[i] = rm_2[i]
					alpha[i] = alpha_2[i]
				end
				
				projTime += toq()

			end # SGD iterations
			#*========= END SGD LOOP ==========
			
			actualNumIter = iter_S
	
			tic()
				
			#*==== ESTIMATE COST AND LINE VIOLATION PROB =====
			cost_soln[iter_S], cost_power_soln[iter_S], anyLineViolProb_soln[iter_S], sumLineViol_soln[iter_S], cost_reserves_soln::Float64, ~, ~, ~, ~, ~, avgLineFlowViol, maxLineFlowViol = estimateCostOfSoln(p0,rp,rm,alpha,tau_sat,gamma_res,tau_pos,gamma_line_iter,delta_line,false,false)
			#*====================================
			
			objTime += toq()

			# update solution for this replicate
			for j = 1:numGen
				p0_solution[iter_S,j] = p0[j]
				rp_solution[iter_S,j] = rp[j]
				rm_solution[iter_S,j] = rm[j]
				alpha_solution[iter_S,j] = alpha[j]
			end

			# update best found solution	
			if(cost_soln[iter_S] < best_cost)
				best_iter = iter_S
				best_cost = cost_soln[iter_S]
				best_cost_power = cost_power_soln[iter_S]
				best_lineViolProb = anyLineViolProb_soln[iter_S]
				best_sumLineViol = sumLineViol_soln[iter_S]
				for j = 1:numGen
					p0_best_solution[j] = p0[j]
					rp_best_solution[j] = rp[j]
					rm_best_solution[j] = rm[j]
					alpha_best_solution[j] = alpha[j]
				end
			end
	
			solutionTime = toq()

			#*========= STORE RESULTS ==========
			if(storeResults)
				
				open(cost_file, "a") do f
					write(f,"$(cost_soln[iter_S]) \n")
				end
				
				open(cost_power_file, "a") do f
					write(f,"$(cost_power_soln[iter_S]) \n")
				end
				
				open(cost_reserves_file, "a") do f
					write(f,"$cost_reserves_soln \n")
				end
				
				open(sum_viol_file, "a") do f
					write(f,"$(sumLineViol_soln[iter_S]) \n")
				end
				
				open(viol_prob_file, "a") do f
					write(f,"$(anyLineViolProb_soln[iter_S]) \n")
				end
				
				open(avg_viol_file, "a") do f
					write(f,"ITER NUM. $iter_S \n")
					for k = 1:numLines
						write(f,"($(lines[k][1]),$(lines[k][2])): $(avgLineFlowViol[k])/$(maxLineFlowViol[k])  \n")
					end
					write(f,"\n")
				end

				open(soln_file, "a") do f
					write(f,"ITER NUM. $iter_S \n")
					write(f,"Expected cost: $(cost_soln[iter_S]) \n")
					write(f,"Expected cost without line penalty costs: $(cost_power_soln[iter_S]) \n")
					write(f,"Expected sum of line flow limit violations: $(sumLineViol_soln[iter_S]) \n")
					write(f,"ANY line limit violation prob: $(anyLineViolProb_soln[iter_S]) \n")
					write(f,"p0: $(p0_solution[iter_S,:]) \n")
					write(f,"rp: $(rp_solution[iter_S,:]) \n")
					write(f,"rm: $(rm_solution[iter_S,:]) \n")
					write(f,"alpha: $(alpha_solution[iter_S,:]) \n")
					write(f,"\n\n")
				end

				open(time_file, "a") do f
					write(f,"$solutionTime \n")
				end
			end
			#*==================================

			if(printIterationStats)
				println("EXPECTED VALUES:  COST (incl. pen.): ",round.(cost_soln[iter_S],displayPrecision),"  COST OF POWER ALONE (incl. reserve pen.): ",round.(cost_power_soln[iter_S],displayPrecision),"  SUM OF LINE VIOL: ",round.(sumLineViol_soln[iter_S],displayPrecision),"  PROB ANY LINE LIMIT VIOLATED: ",round.(anyLineViolProb_soln[iter_S],displayPrecision))
				
				otherTime::Float64 = solutionTime - gradTime - projTime - objTime
				
				println("SOLN TIME: ",round.(solutionTime,2)," GRAD TIME: ",round.(gradTime,2)," PROJ TIME: ",round.(projTime,2)," OBJ TIME: ", round.(objTime,2)," OTHER TIME: ", round.(otherTime,2))
			end

			#*========= HEURISTICS FOR UPDATING STEP LENGTH ==========
			
			# if the min. number of replicates are done, check if termination criteria are satisfied	
			if(iter_S >= minNumReplicates)
				improvements = zeros(Float64,termCheckPeriod)
				for j = 1:termCheckPeriod
					improvements[j] = local_best_cost - cost_soln[iter_S+1-j]
				end
				
				localImpReq::Float64 = impReqForTerm*local_best_cost
				max_improvement::Float64 = maximum(improvements)
				
				if(max_improvement < localImpReq)
					break
				end
			end
			
			# once every stepLengthDecrCheckIter iterations, see if the step length can be decreased	
			if(rem(iter_S,stepLengthDecrCheckIter) == 0)
				improvements = zeros(Float64,stepLengthDecrCheckIter)
				for j = 1:stepLengthDecrCheckIter
					improvements[j] = local_best_cost2 - cost_soln[iter_S+1-j]
				end
				
				localImpReq3::Float64 = impReqForDecrStepSize*local_best_cost2
				
				max_improvement = maximum(improvements)			
				if(max_improvement < -localImpReq3)
					stepLength /= stepLengthDecrFactor
					numStepSizeDecr += 1
				end
			end
			
			# once every stepLengthIncrCheckIter iterations, see if the step length can be increased
			if(rem(iter_S,stepLengthIncrCheckIter) == 0)
				improvements = zeros(Float64,stepLengthIncrCheckIter)
				for j = 1:stepLengthIncrCheckIter
					improvements[j] = local_best_cost2 - cost_soln[iter_S+1-j]
				end
				
				localImpReq2::Float64 = impReqForIncrStepSize*local_best_cost3
				
				max_improvement = maximum(improvements)
				min_improvement = minimum(improvements)
				
				if(min_improvement > -localImpReq2)
					if(max_improvement < localImpReq2)
						stepLength *= stepLengthIncrFactor
						numStepSizeIncr += 1
					end
				end
			end

			# update the best found solutions used for termination and step length modification criteria
			if(iter_S > termCheckPeriod)
				if(cost_soln[iter_S-termCheckPeriod] < local_best_cost)
					local_best_cost = cost_soln[iter_S-termCheckPeriod]
				end
			end
			if(iter_S > stepLengthDecrCheckIter)
				if(cost_soln[iter_S-stepLengthDecrCheckIter] < local_best_cost2)
					local_best_cost2 = cost_soln[iter_S-stepLengthDecrCheckIter]
				end
			end
			if(iter_S > stepLengthIncrCheckIter)
				if(cost_soln[iter_S-stepLengthIncrCheckIter] < local_best_cost3)
					local_best_cost3 = cost_soln[iter_S-stepLengthIncrCheckIter]
				end
			end
			#*========================================================
			
		end # S
		#*========= END OPTIMIZATION LOOP ==========

		println("")
		println("Number of step size incr/decr: ",numStepSizeIncr,"/",numStepSizeDecr)
		println("Step length change factor: ",round.(stepLength./initialStepLength,displayPrecision),"\n")


		# extract the best found solution for the regular and wind generators
		const p0_reg_soln = zeros(Float64,numRegGen)
		const rp_reg_soln = zeros(Float64,numRegGen)
		const rm_reg_soln = zeros(Float64,numRegGen)
		const alpha_reg_soln = zeros(Float64,numRegGen)

		const p0_wind_soln = zeros(Float64, numWindGen)
		const rp_wind_soln = zeros(Float64, numWindGen)
		const rm_wind_soln = zeros(Float64, numWindGen)
		const alpha_wind_soln = zeros(Float64, numWindGen)

		for i = 1:numRegGen
			p0_reg_soln[i] = p0_best_solution[i]
			rp_reg_soln[i] = rp_best_solution[i]
			rm_reg_soln[i] = rm_best_solution[i]
			alpha_reg_soln[i] = alpha_best_solution[i]
		end

		for i = 1:numWindGen
			p0_wind_soln[i] = p0_best_solution[i+numRegGen]
			rp_wind_soln[i] = rp_best_solution[i+numRegGen]
			rm_wind_soln[i] = rm_best_solution[i+numRegGen]
			alpha_wind_soln[i] = alpha_best_solution[i+numRegGen]
		end

		tic()
		
		best_cost, best_cost_power, best_lineViolProb, best_sumLineViol, best_cost_reserves::Float64, best_cost_pos_res_pen::Float64, best_cost_neg_res_pen::Float64, power_gen_soln, wind_spill_soln, ~, ~, ~ = estimateCostOfSoln(p0_best_solution,rp_best_solution,rm_best_solution,alpha_best_solution,tau_sat,gamma_res,tau_pos,gamma_line_iter,delta_line,true,true)
		
		finalCostTime::Float64 = toq()


		#*========= PRINT BEST FOUND SOLUTION ==========


		println("\n*********************")
		println("Lowest cost solution corresponding to iteration #",best_iter)
		println("Expected cost: ",round.(best_cost,displayPrecision))
		println("Expected cost without line penalty costs: ",round.(best_cost_power,displayPrecision))
		println("Expected sum of line flow limit violations: ",round.(best_sumLineViol,displayPrecision))
		println("ANY line limit violation prob: ",round.(best_lineViolProb,displayPrecision))
		
		if(printFinalSolution)
			println("Regular generator variables come first")
			println("p0: ",round.(p0_best_solution,displayPrecision))
			println("rp: ",round.(rp_best_solution,displayPrecision))
			println("rm: ",round.(rm_best_solution,displayPrecision))
			println("alpha: ",round.(alpha_best_solution,displayPrecision),"\n")

			println("\n")
			println("Regular generator variables")
			println("p0_reg: ",p0_reg_soln)
			println("rp_reg: ",rp_reg_soln)
			println("rm_reg: ",rm_reg_soln)
			println("alpha_reg: ",alpha_reg_soln)

			println("\n")
			println("Wind generator variables")
			println("p0_wind: ",p0_wind_soln)
			println("rp_wind: ",rp_wind_soln)
			println("rm_wind: ",rm_wind_soln)
			println("alpha_wind: ",alpha_wind_soln)
		end	
		
		println("")
		println("Final cost time: ",round.(finalCostTime,2))		
		println("")

		iterTime = toq()
		println("RAN ",actualNumIter," ITERATIONS WITH TOTAL TIME FOR LINE PENALTY ITERATION: ",round.(iterTime,2))

		
		#*========= STORE RESULTS ==========
		if(storeResults)

			open(best_soln_file, "a") do f
				write(f,"NOTE: regular generator variables come first \n\n")
				write(f,"Expected cost: $best_cost \n")
				write(f,"Expected cost without line penalty costs: $best_cost_power \n")
				write(f,"Expected sum of line flow limit violations: $best_sumLineViol \n")
				write(f,"Probability that ANY line limit is violated: $best_lineViolProb \n\n")
				write(f,"p0: $p0_best_solution \n")
				write(f,"rp: $rp_best_solution \n")
				write(f,"rm: $rm_best_solution \n")
				write(f,"alpha: $alpha_best_solution \n\n")
				write(f,"Regular generator variables: \n")
				write(f,"p0_reg: $p0_reg_soln \n")
				write(f,"rp_reg: $rp_reg_soln \n")
				write(f,"rm_reg: $rm_reg_soln \n")
				write(f,"alpha_reg: $alpha_reg_soln \n \n")
				write(f,"Wind generator variables: \n")
				write(f,"p0_wind: $p0_wind_soln \n")
				write(f,"rp_wind: $rp_wind_soln \n")
				write(f,"rm_wind: $rm_wind_soln \n")
				write(f,"alpha_wind: $alpha_wind_soln \n")
				write(f,"\n")
			end
			
			open(p0_best_soln_file, "a") do f
				for i = 1:numGen
					write(f,"$(p0_best_solution[i]) ")
				end
				write(f,"\n")
			end
			
			open(rp_best_soln_file, "a") do f
				for i = 1:numGen
					write(f,"$(rp_best_solution[i]) ")
				end
				write(f,"\n")
			end
			
			open(rm_best_soln_file, "a") do f
				for i = 1:numGen
					write(f,"$(rm_best_solution[i]) ")
				end
				write(f,"\n")
			end
			
			open(alpha_best_soln_file, "a") do f
				for i = 1:numGen
					write(f,"$(alpha_best_solution[i]) ")
				end
				write(f,"\n")
			end
				
			open(best_cost_file, "a") do f
				write(f,"$best_cost \n")
			end
			
			open(best_cost_power_file, "a") do f
				write(f,"$best_cost_power \n")
			end
			
			open(best_cost_reserves_file, "a") do f
				write(f,"$best_cost_reserves \n")
			end
			
			open(best_cost_pos_res_file, "a") do f
				write(f,"$best_cost_pos_res_pen \n")
			end
			
			open(best_cost_neg_res_file, "a") do f
				write(f,"$best_cost_neg_res_pen \n")
			end
			
			open(best_sum_viol_file, "a") do f
				write(f,"$best_sumLineViol \n")
			end
			
			open(best_viol_prob_file, "a") do f
				write(f,"$best_lineViolProb \n")
			end
			
			open(power_gen_soln_file, "a") do f
				for i = 1:numGen
					write(f,"$(power_gen_soln[i]) ")
				end
				write(f,"\n")
			end
			
			open(wind_spill_soln_file, "a") do f
				for i = 1:numWindGen
					write(f,"$(wind_spill_soln[i]) ")
				end
				write(f,"\n")
			end

			open(time_file, "a") do f
				write(f,"\n\n")
				write(f,"TOTAL TIME: $iterTime")
			end
		end
		#*==================================

		# set initial guess for next iteration as best solution of this iteration
		for i = 1:numGen
			p0_init[i] = p0_best_solution[i]
			rp_init[i] = rp_best_solution[i]
			rm_init[i] = rm_best_solution[i]
			alpha_init[i] = alpha_best_solution[i]
		end

	end # line penalties

	trialTime::Float64 = toq()
	println("\nTime for trial #",trial,": ",round.(trialTime,2))

end # trials

toc() # overall time