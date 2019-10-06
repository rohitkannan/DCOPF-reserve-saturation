# ROUTINES FOR COMPUTING STOCHASTIC GRADIENT OF THE OBJECTIVE FUNCTION


# get the partial derivatives of recourse variables wrt the first-stage variables for a given scenario
function getSensitivities(alpha::Array{Float64},pt::Array{Float64},omega_d::Array{Float64},tau_sat::Float64)

# store values of the derivative of the saturation function with respect to its argument							
	dsat_dpt = zeros(Float64,numGen)
	for i = 1:numGen
		dsat_dpt[i] = smoothSaturationDeriv(pt[i],pmin[i],pmax[i],tau_sat)
	end
	den_s::Float64 = sum(alpha.*dsat_dpt)
		
# first, get partial derivatives wrt the p0s	
	deriv_s_p0 = -dsat_dpt/den_s
	dp_dp0 = (dsat_dpt.*alpha)*(deriv_s_p0') + diagm(dsat_dpt)

	return dp_dp0
end



# compute a stochastic gradient of the objective function wrt the first-stage variables
function getStochasticGradient(p0::Array{Float64},rp::Array{Float64},rm::Array{Float64},alpha::Array{Float64},
								numSamples::Int64,tau_sat::Float64,gamma_res::Float64,tau_pos::Float64,
								gamma_line::Float64,delta_line::Float64)

# generate scenarios for computing the stochastic gradient
	omega_d, omega_w = generateScenarios(numSamples)

	stochgrad_p0 = zeros(Float64,numGen)
	stochgrad_rp = zeros(Float64,numGen)
	stochgrad_rm = zeros(Float64,numGen)
	stochgrad_alpha = zeros(Float64,numGen)
	
	for iter = 1:numSamples
	# first, solve the recourse system. this routine also updates the wind generator levels at the current scenario
		pt, p, slack, theta = solveRecourseProblem(p0,alpha,omega_d[:,iter],omega_w[:,iter],true,tau_sat)
		
	# then, compute partial derivatives of recourse variables wrt first-stage variables at above solution
	# ignore the ones we know to be zero
		dp_dp0 = getSensitivities(alpha,pt,omega_d,tau_sat)
		
	# compute stochastic gradient contribution from this scenario
	# add on the derivatives from the penalties for positive, negative reserves for generator i
		derivPosViol = zeros(Float64,numGen)
		derivNegViol = zeros(Float64,numGen)
		for i = 1:numGen
			posResViol::Float64 = p[i]-p0[i]-rp[i]
			derivPosViol[i] = reserveSmoothPenDeriv(posResViol,posResCost[i],gamma_res,tau_pos)
			
			negResViol::Float64 = p0[i]-p[i]-rm[i]
			derivNegViol[i] = reserveSmoothPenDeriv(negResViol,negResCost[i],gamma_res,tau_pos)
		end
		
		sumDemandDeviations::Float64 = sum(omega_d[j,iter] for j = 1:numLoads) - sum(meanDemand[j] for j = 1:numLoads)
		stochgrad_p0_iter1 = (dp_dp0')*(derivPosViol-derivNegViol)
		stochgrad_p0 += stochgrad_p0_iter1 + (derivNegViol - derivPosViol)
		stochgrad_alpha += stochgrad_p0_iter1*(sumDemandDeviations + slack)
		stochgrad_rp -= derivPosViol
		stochgrad_rm -= derivNegViol

	# add on the derivatives from the line limit penalties for line k
		derivLineViol = zeros(Float64,numLines)
		for k = 1:numLines
			lineFlow::Float64 = beta[k]*(theta[lines[k][1]] - theta[lines[k][2]])
			derivLineViol[k] = linePenDeriv(lineFlow,flow_max[k],gamma_line,delta_line)
		end	
		
		stochgrad_p0_iter2 = (dp_dp0')*(Line_mat_times_B_over_G_trans*(derivLineViol.*beta))
		stochgrad_p0 += stochgrad_p0_iter2
		stochgrad_alpha += stochgrad_p0_iter2*(sumDemandDeviations + slack)
	end
	
	stochgrad_p0 /= numSamples
	stochgrad_rp /= numSamples
	stochgrad_rm /= numSamples
	stochgrad_alpha /= numSamples

# add the derivatives from the first-stage generation costs
	stochgrad_p0 += genCost
	stochgrad_rp += posResCost
	stochgrad_rm += negResCost
	
	return stochgrad_p0, stochgrad_rp, stochgrad_rm, stochgrad_alpha
end