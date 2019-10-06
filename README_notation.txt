# SUMMARY OF THE SOME KEY NOTATION

# numBuses: number of buses
# numRegGen: number of regular generators
# numWindGen: number of wind generators
# numGen: total number of regular + wind generators
# numLoads: number of load/demand nodes

# regGenToBus: indicate bus number that a regular generator serves
# windGenToBus: indicate bus number that a wind generator serves
# loadToBus: indicate bus number at which the load is required

# regGenAtBus: which regular generators cater to the bus?
# windGenAtBus: which wind generators cater to the bus?
# loadAtBus: which loads need to be met at the bus?

# pmin: minimum power outputs of the regular + wind generators
# pmax: maximum power outputs of the regular + wind generators
# pmax_det: maximum power outputs of the regular + wind generators using deterministic wind predictions
# pmax_FS: maximum power outputs of the regular + wind generators that are used to make first-stage decisions

# genCost: cost of power from the regular + wind generators
# posResCost: cost of up/positive reserves for the regular + wind generators
# negResCost: cost of down/negative reserves for the regular + wind generators

# genRes: which regular + wind generators are ready to provide reserves?
# genGuarRes: which regular + wind generators are ALWAYS READY to provide reserves?
# min_alpha: minimum participation factors for the generators in genGuarRes

# lines: lines between buses in the network
# beta: line susceptances
# flow_max: line flow limits
# B_mat: matrix multiplying the theta variables in the DCOPF equation system

# omega_d: load/demand scenarios
# omega_w: wind output scenarios
# meanDemand: average load demands
# demandStdev: standard deviation of load demands
# windPredMean: average wind power output predictions
# windPredStdev: standard deviation of wind power output predictions

# gamma_line: line flow limit violation penalty parameter
# delta_line: beyond what fraction of the line flow limit do we start penalizing line flows?
# gamma_res: penalty factor for exceeding up/down reserves
# gamma_gen: penalty factor for generation levels exceeding generator limits
# genBoundEps: chance constraint violation probability allowance for the CAP model
# tau_sat: relative smoothing parameter for the smooth saturation function
# tau_pos: smoothing parameter for the smooth positive function
# slackTolerance: tolerance for the bisection procedure used to solve the recourse system
# windBoundFactor: factor for determining pmax_FS for the wind generators


# NOTE: the solution for the regular generators are stored first, followed by the solution for the wind generators

# p0: first-stage power generation decisions for the regular + wind generators
# rp: first-stage up/positive reserve decisions for the regular + wind generators
# rm: first-stage down/negative reserve decisions for the regular + wind generators
# alpha: first-stage participation factor decisions for the regular + wind generators

# pt: recourse power target decisions for the regular + wind generators
# p: recourse power generation decisions for the regular + wind generators
# theta: recourse bus phase angle decisions
# slack: slack variable that accounts for generator power output saturation while using affine control policy