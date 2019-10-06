# ROUTINES FOR EVALUATING QUALITY OF COMPUTED FIRST-STAGE SOLUTION


#*====== CONSTRUCT SCENARIOS FOR MC ESTIMATION =======
# first, generate scenarios for estimating progress of stochastic approximation method
const numMCScenarios = 10000

srand(7613)
omega_d_MC, omega_w_MC = generateScenarios(numMCScenarios)
srand()


# then, generate an independent batch of scenarios for evaluating quality of solutions of ALL METHODS
const numAnalyticalScenarios = 100000

srand(5338)
omega_d_analytical, omega_w_analytical = generateScenarios(numAnalyticalScenarios)
srand()



# function to estimate (expected) cost of first-stage decisions
function estimateCostOfSoln(p0::Array{Float64},rp::Array{Float64},rm::Array{Float64},alpha::Array{Float64},
								tau_sat::Float64,gamma_res::Float64,tau_pos::Float64,gamma_line::Float64,delta_line::Float64,
								useAnalyticalScen::Bool=false,printSoln::Bool=false,printLineViolStats::Bool=false)
							
# expected first-stage + recourse costs (including all penalties)
	expectedCost::Float64 = 0.0
# expected cost of power purchasing decisions (does not include line penalties, but includes reserve penalties)
	expectedCostOfPower::Float64 = 0.0
# expected cost of reserves (including reserve penalties)
	expectedCostOfReserves::Float64 = 0.0
# expected positive reserve penalty
	expectedPosResPenalty::Float64 = 0.0
# expected negative reserve penalty
	expectedNegResPenalty::Float64 = 0.0
	
# estimate of probability that ANY line flow limit is violated
	probAnyLineFlowViol::Float64 = 0.0
# expected sum of line flow limit violations
	sumLineFlowViol::Float64 = 0.0
	
# estimate of probability that each line flow limit is violated
	lineFlowViolProb = zeros(Float64,numLines)
# average violation of each line flow limit (averaged over the violated cases)
	avgAbsLineFlowViol = zeros(Float64,numLines)
# maximum violation of each line flow limit
	maxAbsLineFlowViol = zeros(Float64,numLines)
	
# expected amount of wind power generation
	expectedPowerGen = zeros(Float64,numGen)
# expected amount of wind power not utilized
	expectedWindSpill = zeros(Float64,numWindGen)
	
	numScenarios::Int64 = numMCScenarios
	if(useAnalyticalScen)
		numScenarios = numAnalyticalScenarios
	end

# loop over scenarios for MC estimation	
	for scen = 1:numScenarios
		p = zeros(Float64,numGen)
		
		if(useAnalyticalScen)
			~, p, ~, theta = solveRecourseProblem(p0,alpha,omega_d_analytical[:,scen],omega_w_analytical[:,scen],false,tau_sat)
		else
			~, p, ~, theta = solveRecourseProblem(p0,alpha,omega_d_MC[:,scen],omega_w_MC[:,scen],false,tau_sat)
		end
		
	# compute expected power generation quantities
		for i = 1:numGen
			expectedPowerGen[i] += p[i]
		end

	# calculate expected wind spill
		for i = 1:numWindGen
			if(useAnalyticalScen)
				expectedWindSpill[i] += (omega_w_analytical[i,scen] - p[i+numRegGen])
			else
				expectedWindSpill[i] += (omega_w_MC[i,scen] - p[i+numRegGen])
			end
		end
		
	# add reserve penalty recourse costs for positive and negative reserves
		for i = 1:numGen
			posResViol::Float64 = p[i]-p0[i]-rp[i]
			negResViol::Float64 = p0[i]-p[i]-rm[i]
		
			posResPenCost::Float64 = reserveSmoothPen(posResViol,posResCost[i],gamma_res,tau_pos)
			negResPenCost::Float64 = reserveSmoothPen(negResViol,negResCost[i],gamma_res,tau_pos)
			
			expectedCost += posResPenCost + negResPenCost
			expectedCostOfPower += posResPenCost + negResPenCost
			expectedCostOfReserves += posResPenCost + negResPenCost
			expectedPosResPenalty += posResPenCost
			expectedNegResPenalty += negResPenCost
		end
	
	# add line penalty costs and calculate probability of line flow violation
		anyLineViol::Bool = false
		for k = 1:numLines
			lineFlow::Float64 = beta[k]*(theta[lines[k][1]] - theta[lines[k][2]])
			
			expectedCost += linePen(lineFlow,flow_max[k],gamma_line,delta_line)
			
			if(abs(lineFlow) > flow_max[k])
				lineFlowViolProb[k] += 1.0
				anyLineViol = true
				
				lineFlowViol::Float64 = (abs(lineFlow) - flow_max[k])
				sumLineFlowViol += lineFlowViol
				avgAbsLineFlowViol[k] += lineFlowViol
				if(lineFlowViol > maxAbsLineFlowViol[k])
					maxAbsLineFlowViol[k] = lineFlowViol
				end
			end
		end
		
		if(anyLineViol)
			probAnyLineFlowViol += 1.0
		end
	end


# re-scale the average line flow violations and Monte Carlo estimates	
	for k = 1:numLines
		if(lineFlowViolProb[k] > 0)
			avgAbsLineFlowViol[k] /= lineFlowViolProb[k]
		end
	end
	
	expectedCost /= numScenarios
	expectedCostOfPower /= numScenarios
	expectedCostOfReserves /= numScenarios
	expectedPosResPenalty /= numScenarios
	expectedNegResPenalty /= numScenarios
	probAnyLineFlowViol /= numScenarios
	lineFlowViolProb /= numScenarios
	sumLineFlowViol /= numScenarios
	expectedPowerGen /= numScenarios
	expectedWindSpill /= numScenarios

	
# add first-stage costs	
	for i = 1:numGen
		genPowerCost::Float64 = genCost[i]*p0[i] + posResCost[i]*rp[i] + negResCost[i]*rm[i]
		expectedCost += genPowerCost
		expectedCostOfPower += genPowerCost
		expectedCostOfReserves += (genPowerCost - genCost[i]*p0[i])
	end

	
	if(printSoln)
		println("Expected cost (incl. all penalties): ",round.(expectedCost,displayPrecision))
		println("Expected cost without line penalty costs: ",round.(expectedCostOfPower,displayPrecision))
		println("Expected cost of reserves (incl. reserve penalties): ",round.(expectedCostOfReserves,displayPrecision))
		println("Expected reserve penalty costs: ",round.(expectedPosResPenalty+expectedNegResPenalty,displayPrecision))
		
		println("Expected sum of line flow violations: ",round.(sumLineFlowViol,displayPrecision))
		println("Probability that ANY line flow is violated: ",round.(probAnyLineFlowViol,displayPrecision))
		
	# do we want to print out detailed violation statistics for each line?
	# can be too much output for large instances
		if(printLineViolStats)		
			println("\nProbability that each line flow is violated: ")
			for k = 1:numLines
				print("(",lines[k][1],",",lines[k][2],"): ",round.(lineFlowViolProb[k],displayPrecision),"   ")
			end
			
			println("\n\nAverage/Maximum line flow violations: ")
			for k = 1:numLines
				print("(",lines[k][1],",",lines[k][2],"): ",round.(avgAbsLineFlowViol[k],displayPrecision),"/",
													round.(maxAbsLineFlowViol[k],displayPrecision),"   ")
			end
		end
		println("\n")
	end
	
	return expectedCost, expectedCostOfPower, probAnyLineFlowViol, sumLineFlowViol, expectedCostOfReserves, expectedPosResPenalty, expectedNegResPenalty, expectedPowerGen, expectedWindSpill, lineFlowViolProb, avgAbsLineFlowViol, maxAbsLineFlowViol
end


# estimate the probability of generator bounds violation and other stats
function estimateGenBoundViolStats(p0::Array{Float64},alpha::Array{Float64},printViolProb::Bool=false,
										printAllStats::Bool=false)

# probability that ANY generator bounds are violated
	probAnyGenBoundViol::Float64 = 0.0

# probability that individual generator bounds are violated
	probGenBoundViol = zeros(Float64,numGen)

# maximum generator bounds violation
	maxGenBoundViol = zeros(Float64,numGen)

# median generator bounds violation
	medianGenBoundViol = zeros(Float64,numGen)

# store all the generator bounds violations for all generators and for all scenarios
	genBoundViol = zeros(Float64,numGen,numAnalyticalScenarios)
	
	sumMeanDemands = sum(meanDemand[j] for j=1:numLoads)
	sumDemands = (sum(omega_d_analytical,1))'
	sumDemandDeviations = sumDemands - sumMeanDemands*ones(Float64,numAnalyticalScenarios)
	anyGenBoundViolated::Bool = false
	boundViolTol::Float64 = 1E-09
	powerOutput::Float64 = 0.0
	
	for scen = 1:numAnalyticalScenarios
		# update generator bounds based on the wind scenarios
		updateGenBounds(omega_w_analytical[:,scen])
		
		anyGenBoundViolated = false
		for i = 1:numGen
			powerOutput = p0[i] + alpha[i]*sumDemandDeviations[scen]
			genBoundViol[i,scen] = max(powerOutput - pmax[i], pmin[i] - powerOutput, 0.0)
			if(genBoundViol[i,scen] > boundViolTol)
				anyGenBoundViolated = true
				probGenBoundViol[i] += 1.0
				if(genBoundViol[i,scen] > maxGenBoundViol[i])
					maxGenBoundViol[i] = genBoundViol[i,scen]
				end
			end
		end
		if(anyGenBoundViolated)
			probAnyGenBoundViol += 1.0
		end
	end

	# re-scale probabilities
	probAnyGenBoundViol /= numAnalyticalScenarios
	probGenBoundViol /= numAnalyticalScenarios
	medianGenBoundViol = median(genBoundViol,2)

	
	if(printViolProb)
		println("Probability that ANY generator's bounds are violated by the affine policy: ",round.(probAnyGenBoundViol,displayPrecision))
		
	# do we want to print out detailed violation statistics for each generator?
	# can be too much output for large instances
		if(printAllStats)		
			println("Probability that each generator bound is violated: ",round.(probGenBoundViol,displayPrecision))
			
			println("Max. generator bounds violation: ",round.(maxGenBoundViol,displayPrecision))
			
			println("Median generator bounds violation: ",round.(medianGenBoundViol,displayPrecision))
		end
		println("\n")
	end

	return probAnyGenBoundViol, probGenBoundViol, maxGenBoundViol, medianGenBoundViol
end