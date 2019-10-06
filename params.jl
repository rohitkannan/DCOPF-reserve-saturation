# SET MODEL PARAMETERS


#*========== MAIN PARAMETERS ===========
# beyond what fraction of line limit do we penalize?
const delta_line = 0.95

# how much costlier are the regular generator reserves compared to normal costs?
const reservesCostFactor = 1.5

# penalty parameter for exceeding reserves in the recourse model
const gamma_res = 10.0

# how do the wind generator reserve costs scale compared to the cheapest regular generator?
const windReservesCostFactor = 0.1


#*========== OTHER PARAMETERS ===========
# parameter for smoothing saturation function (relative to bounds)
const tau_sat = 1E-04

# parameter for smooth positive function approximation
const tau_pos = 1E-04

# tolerance for solving recourse system of nonlinear equations
const slackTolerance = 1E-06

# display precision
const displayPrecision = 4

# factor by which to scale wind forecast standard deviation to compute first-stage wind generator upper bounds (pmax_FS)
const windBoundFactor = 5.0