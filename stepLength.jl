# ROUTINES FOR ESTIMATING THE STEP LENGTH FOR THE STOCHASTIC APPROXIMATION ALGORITHM


# get the stochastic gradient of the first-stage objective function wrt the first-stage variables for GIVEN SCENARIOS
function getStochasticGradientForFixedSample(p0::Array{Float64},rp::Array{Float64},rm::Array{Float64},alpha::Array{Float64},
						omega_d::SubArray{Float64,2,Array{Float64,2},Tuple{Base.Slice{Base.OneTo{Int64}},UnitRange{Int64}},true},
						omega_w::SubArray{Float64,2,Array{Float64,2},Tuple{Base.Slice{Base.OneTo{Int64}},UnitRange{Int64}},true},
						tau_sat::Float64,gamma_res::Float64,tau_pos::Float64,gamma_line::Float64,delta_line::Float64)

	numSamples::Int64 = size(omega_d,2)

	stochgrad_p0 = zeros(Float64,numGen)
	stochgrad_rp = zeros(Float64,numGen)
	stochgrad_rm = zeros(Float64,numGen)
	stochgrad_alpha = zeros(Float64,numGen)
	
	for iter = 1:numSamples
	# first, solve the recourse system
		pt, p, slack, theta = solveRecourseProblem(p0,alpha,omega_d[:,iter],omega_w[:,iter],true,tau_sat)

	# then, compute partial derivatives of recourse variables wrt first-stage variables at above solution. ignore the ones we know to be zero
		dp_dp0 = getSensitivities(alpha,pt,omega_d[:,iter],tau_sat)

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



# estimate Lipschitz constant of the gradient of the objective function
function estimateLipschitzConstant(p0_ref::Array{Float64},rp_ref::Array{Float64},rm_ref::Array{Float64},alpha_ref::Array{Float64},
									numSamples::Int64,numGradientSamples::Int64,tau_sat::Float64,gamma_res::Float64,
									tau_pos::Float64,gamma_line::Float64,delta_line::Float64)

# generate scenarios of the random variables								
	omega_d, omega_w = generateScenarios(numGradientSamples*numSamples)

# four sets of Lipschitz constant estimates for the four types of first-stage variables
	L_max = zeros(Float64,4)
	for iter = 1:numSamples	
	# sample points in the first-stage feasible region that are close to the reference point
		p0_1_int, rp_1_int, rm_1_int, alpha_1_int = pickPointOnSphere(p0_ref,rp_ref,rm_ref,alpha_ref)
		p0_1, rp_1, rm_1, alpha_1, solve_status_1 = projectOntoFeasibleSet(p0_1_int,rp_1_int,rm_1_int,alpha_1_int)
		
		p0_2_int, rp_2_int, rm_2_int, alpha_2_int = pickPointOnSphere(p0_ref,rp_ref,rm_ref,alpha_ref)
		p0_2, rp_2, rm_2, alpha_2, solve_status_2 = projectOntoFeasibleSet(p0_2_int,rp_2_int,rm_2_int,alpha_2_int)
		
		if((solve_status_1 == :Optimal || solve_status_1 == :Suboptimal) && 
			(solve_status_2 == :Optimal || solve_status_2 == :Suboptimal))	
			# determine which samples of the random variables will be used to approximate the gradient of the objective function
			omega_start::Int64 = (iter-1)*numGradientSamples + 1
			omega_end::Int64 = iter*numGradientSamples
		
			# approximate gradient of the objective function using stochastic gradients
			@views grad_p0_1, grad_rp_1, grad_rm_1, grad_alpha_1 = getStochasticGradientForFixedSample(p0_1,rp_1,rm_1,alpha_1,
																	omega_d[:,omega_start:omega_end],omega_w[:,omega_start:omega_end],
																	tau_sat,gamma_res,tau_pos,gamma_line,delta_line)
			@views grad_p0_2, grad_rp_2, grad_rm_2, grad_alpha_2 = getStochasticGradientForFixedSample(p0_2,rp_2,rm_2,alpha_2,
																	omega_d[:,omega_start:omega_end],omega_w[:,omega_start:omega_end],
																	tau_sat,gamma_res,tau_pos,gamma_line,delta_line)

			L_est = zeros(Float64,4)
			L_est[1] = norm(grad_p0_2 - grad_p0_1)/norm(p0_2 - p0_1)
			L_est[2] = norm(grad_rp_2 - grad_rp_1)/norm(rp_2 - rp_1)
			L_est[3] = norm(grad_rm_2 - grad_rm_1)/norm(rm_2 - rm_1)
			L_est[4] = norm(grad_alpha_2 - grad_alpha_1)/norm(alpha_2 - alpha_1)
			
			for i = 1:4
				if(L_est[i] > L_max[i])
					L_max[i] = L_est[i]
				end
			end
		end
	end
	
	return L_max
end


# estimate step length for the stochastic approximation method using the simple rule in Ghadimi, Lan, and Zhang
# we estimate a different step length for each of the four sets of variables p0, rp, rm, and alpha because the problem is ill-conditioned
function estimateStepLength(p0_ref::Array{Float64},rp_ref::Array{Float64},rm_ref::Array{Float64},alpha_ref::Array{Float64},
							numGradientSamples::Int64,numSamplesForLipEst::Int64,tau_sat::Float64,gamma_res::Float64,
							tau_pos::Float64,gamma_line::Float64,delta_line::Float64)

	rho_est = estimateLipschitzConstant(p0_ref,rp_ref,rm_ref,alpha_ref,numSamplesForLipEst,numGradientSamples,tau_sat,
													gamma_res,tau_pos,gamma_line,delta_line)
# step length is an array of size four
	stepLength = 1.0./(2.0*rho_est)
	
	return stepLength
end


# given a reference point, pick a random point within a given radius
function pickPointOnSphere(p0_orig::Array{Float64},rp_orig::Array{Float64},rm_orig::Array{Float64},alpha_orig::Array{Float64})

# factor to determine sampling neighborhood
	scaleFactor::Float64 = 10.0

# generate point on the unit sphere in the dimension of the first-stage variables
	unifRand = rand(1)
	factor1::Float64 = (unifRand[1])^(1.0/numGen)
	normRand = rand(Normal(0,1), numGen)
	factor2::Float64 = norm(normRand)
# this is a vector
	factor3 = factor1*normRand/factor2

	p0_samp = p0_orig + (abs.(p0_orig)).*factor3/scaleFactor
	rp_samp = rp_orig + (abs.(rp_orig)).*factor3/scaleFactor
	rm_samp = rm_orig + (abs.(rm_orig)).*factor3/scaleFactor
	alpha_samp = alpha_orig + (abs.(alpha_orig)).*factor3/scaleFactor

	return p0_samp, rp_samp, rm_samp, alpha_samp
end