# DEFINE FUNCTIONS THAT ARE PART OF THE MODEL



#*====== SMOOTH SATURATION FUNCTION =======

# true nonsmooth saturation function
function saturation(x::Float64,x_min::Float64,x_max::Float64)

	x_out::Float64 = median([x, x_min, x_max])
	return x_out
end


# smooth saturation function
function smoothSaturation(x::Float64,x_min::Float64,x_max::Float64,tau_sat::Float64)

# scale smoothing parameter relative to x_max and x_min	
	tau::Float64 = tau_sat*(x_max - x_min)

	x_out::Float64 = x_max

	if(x < x_min - tau)
		x_out = x_min
	elseif(x <= x_min + tau)
		x_out = x_min + (1/(4*tau))*(x - (x_min - tau))^2
	elseif(x < x_max - tau)
		x_out = x
	elseif(x <= x_max + tau)
		x_out = x_max - (1/(4*tau))*(x - (x_max + tau))^2
	end

	return x_out
end


# return derivative of smooth saturation function
function smoothSaturationDeriv(x::Float64,x_min::Float64,x_max::Float64,tau_sat::Float64)

# scale smoothing parameter relative to x_max and x_min	
	tau::Float64 = tau_sat*(x_max - x_min)
	
	deriv::Float64 = 0

	if(x < x_min - tau)
		# do nothing
	elseif(x <= x_min + tau)
		deriv = (1/(2*tau))*(x - (x_min - tau))
	elseif(x < x_max - tau)
		deriv = 1
	elseif(x <= x_max + tau)
		deriv = - (1/(2*tau))*(x - (x_max + tau))
	end

	return deriv
end


#*====== LINE PENALTY FUNCTION =======

# line penalty function
function linePen(x::Float64,x_max::Float64,gamma_line::Float64,delta_line::Float64)

	line_viol::Float64 = maximum([abs(x) - delta_line*x_max, 0])
	pen_val::Float64 = gamma_line*(line_viol)^2
	return pen_val
end


# derivative of line penalty function
function linePenDeriv(x::Float64,x_max::Float64,gamma_line::Float64,delta_line::Float64)

	line_viol::Float64 = maximum([abs(x) - delta_line*x_max, 0])
	deriv::Float64 = 2*gamma_line*line_viol*sign(x)
	return deriv
end


#*====== RESERVE PENALTY FUNCTION =======

# returns smooth version of max{x, 0}
function smoothPos(x::Float64,tau_pos::Float64)

	if(x <= 0)
		return (tau_pos*log(1+exp(x/tau_pos)))
	else
	# rewrite for numerical stability
		return (x + tau_pos*log(1 + exp(-x/tau_pos)))
	end

end


# derivative of smooth positive function
function smoothPosDeriv(x::Float64,tau_pos::Float64)

	deriv::Float64 = 1/(1+exp(-x/tau_pos))
	return deriv
end


# reserves penalty function
function reserveSmoothPen(x::Float64,costCoeff::Float64,gamma_res::Float64,tau_pos::Float64)
	
	res_pen::Float64 = gamma_res*costCoeff*smoothPos(x,tau_pos)
	return res_pen
end


# derivative of reserves penalty function
function reserveSmoothPenDeriv(x::Float64,costCoeff::Float64,gamma_res::Float64,tau_pos::Float64)
	
	deriv::Float64 = gamma_res*costCoeff*smoothPosDeriv(x,tau_pos)
	return deriv
end