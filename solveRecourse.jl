# SOLVE THE RECOURSE PROBLEM FOR GIVEN VALUES OF THE FIRST-STAGE VARIABLES AND A GIVEN SCENARIO


# solve recourse problem by bisection
function solveRecourseProblem(p0::Array{Float64},alpha::Array{Float64},omega_d::Array{Float64},
									omega_w::Array{Float64},useSmoothSatFn::Bool,tau_sat::Float64)

# update generator bounds based on the wind scenarios
	updateGenBounds(omega_w)

# compute demand deviations for control policy	
	sumDemands::Float64 = sum(omega_d[j] for j = 1:numLoads)
	sumDemandDeviations::Float64 = sumDemands - sum(meanDemand[j] for j = 1:numLoads)
	
# check if the recourse problem has a solution
	demandBalance_lower::Float64 = -sumDemands
	demandBalance_upper::Float64 = -sumDemands
	for i = 1:numGen
		if(alpha[i] > 0.0)
			demandBalance_lower += pmin[i]
			demandBalance_upper += pmax[i]
		else
			demandBalance_lower += median([p0[i],pmin[i],pmax[i]])
			demandBalance_upper += median([p0[i],pmin[i],pmax[i]])
		end
	end
	
	if(demandBalance_upper < 0.0 || demandBalance_lower > 0.0)
		println("p0: ",p0,"\nalpha: ",alpha,"\npmin: ",pmin,"\npmax: ",pmax,"\nsumDemands: ",sumDemands)
		throw(ErrorException("Relatively complete recourse does not hold!!!"))
	end

	pt = p0 + alpha*sumDemandDeviations
	p = zeros(Float64,numGen)

# check if the affine policy is adequate, or if we need to provide up/down reserves	
# also compute the power generation quantities p using the saturation function
	for i = 1:numGen
		if(useSmoothSatFn)
			p[i] = smoothSaturation(pt[i],pmin[i],pmax[i],tau_sat)
		else
			p[i] = saturation(pt[i],pmin[i],pmax[i])
		end
	end
	netLoadDeviation::Float64 = sum(p) - sumDemands
	
# minimum participation factor to determine bounds on the slack
	minPartFactor::Float64 = min_alpha/11.0

# determine if we need to provide up/down reserves, and get trivial bounds on slack variable
# if netLoadDeviation < 0, we need up reserves, if > 0, we need down reserves
	slack_lower::Float64 = 0.0
	slack_upper::Float64 = 0.0	
	tmpVar::Float64 = 0.0
	tau_1::Float64 = 0.0
	if(netLoadDeviation < 0.0)
		slack_lower = -netLoadDeviation
		slack_upper = -netLoadDeviation
		for i = 1:numGen
			if(p[i] >= pmax[i])
			# do nothing
			elseif(alpha[i] >= minPartFactor)
				if(useSmoothSatFn)
					tau_1 = tau_sat*(pmax[i] - pmin[i])
					tmpVar = (pmax[i] + tau_1 - p0[i] - alpha[i]*sumDemandDeviations)/alpha[i]
				else
					tmpVar = (pmax[i] - p0[i] - alpha[i]*sumDemandDeviations)/alpha[i]
				end
				if(tmpVar > slack_upper)
					slack_upper = tmpVar
				end
			end
		end
	elseif(netLoadDeviation > 0.0)
		slack_upper = -netLoadDeviation
		slack_lower = -netLoadDeviation
		for i = 1:numGen
			if(p[i] <= pmin[i])
			# do nothing
			elseif(alpha[i] >= minPartFactor)
				if(useSmoothSatFn)
					tau_1 = tau_sat*(pmax[i] - pmin[i])
					tmpVar = (pmin[i] - tau_1 - p0[i] - alpha[i]*sumDemandDeviations)/alpha[i]
				else
					tmpVar = (pmin[i] - p0[i] - alpha[i]*sumDemandDeviations)/alpha[i]
				end
				if(tmpVar < slack_lower)
					slack_lower = tmpVar
				end
			end
		end
	end

	
	slack_val::Float64 = (slack_lower + slack_upper)/2.0
	absoluteDemandError::Float64 = 100.0
	relativeDemandError::Float64 = 100.0
	netLoadDeviation_iter::Float64 = 0.0
	sum_p::Float64 = 0.0
	
# solve recourse problem by bisection by accounting for reserves
	numIter::Int64 = 0
	while(true)
		numIter += 1
		pt = p0 + alpha*(sumDemandDeviations + slack_val)
		for i = 1:numGen
			if(useSmoothSatFn)
				p[i] = smoothSaturation(pt[i],pmin[i],pmax[i],tau_sat)
			else
				p[i] = saturation(pt[i],pmin[i],pmax[i])
			end
		end
		
		sum_p = sum(p)
		netLoadDeviation_iter = sum_p - sumDemands
		absoluteDemandError = abs(netLoadDeviation_iter)
		relativeDemandError = absoluteDemandError/max(abs(sum_p),abs(sumDemands))
		
		if(absoluteDemandError <= slackTolerance || relativeDemandError <= slackTolerance)
			break
		end

		if(netLoadDeviation_iter > 0.0)
			slack_upper = slack_val - netLoadDeviation_iter
		elseif(netLoadDeviation_iter < 0.0)
			slack_lower = slack_val - netLoadDeviation_iter
		end	
		
		if(numIter > 1000)
			println("***  ERROR IN solveRecourseProblem  ***")
			println("useSmoothSatFn: ",useSmoothSatFn,"\np0: ",p0,"\nalpha: ",alpha,"\npmin: ",pmin,"\npmax: ",pmax,"\nsumDemands: ",sumDemands,"\nsumDemandDeviations: ",sumDemandDeviations,"\nnetLoadDeviation: ",netLoadDeviation)
			println("slack_lower: ",slack_lower,"  slack_upper: ",slack_upper,"  slack_val: ",slack_val,"  netLoadDeviation_iter: ",netLoadDeviation_iter,"  absoluteDemandError: ",absoluteDemandError,"  relativeDemandError: ",relativeDemandError)
			throw(ErrorException("Too many bisection iterations in solveRecourseProblem!"))
		end
		
		slack_val = (slack_lower + slack_upper)/2.0
	end
	
	sum_p = sum(p)
	relativeDemandError = abs(sum_p - sumDemands)/max(abs(sumDemands),abs(sum_p))
	absoluteDemandError = abs(sum_p - sumDemands)
	
# final check as to whether the solution for p is sufficiently accurate	
	if(relativeDemandError > 10.0*slackTolerance && absoluteDemandError > 10.0*slackTolerance)
		println("***  ERROR IN solveRecourseProblem  ***")
		println("sum generated power: ",sum_p)
		println("sumDemands: ",sumDemands)
		println("absoluteDemandError: ",absoluteDemandError,"  relativeDemandError: ",relativeDemandError,"  tolerance: ", slackTolerance)
		throw(ErrorException("Recourse system was not solved to enough accuracy! Check that relatively complete recourse holds"))
	end
	
	
# finally, solve system of DCOPF equations for the phase angles theta
	theta = zeros(Float64,numBuses)
	theta[1:numBuses-1] = B_over_G*p - B_over_L*omega_d

	return pt, p, slack_val, theta
end