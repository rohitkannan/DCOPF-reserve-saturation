# 6 Bus Model based on http://motor.ece.iit.edu/data/6bus_Data_ES.pdf

# OPTIONS for caseID: 
# wind_reserves: wind generator can provide reserves. this is the default
# no_wind_reserves: wind generator cannot provide reserves
# nondispatchable_wind: wind generators are not dispatchable

const caseID = "no_wind_reserves"


# include basic model parameters
include("params.jl")



#*========== BASIC TOPOLOGY ===========

# number of buses, generators (regular and wind), loads
const numBuses = 6
const numGen = 3
const numRegGen = 2
const numWindGen = 1
const numLoads = 3


# map regular generators to buses
const regGenToBus = [1,6]

# map wind generators to buses
const windGenToBus = [2]

# map loads to buses
const loadToBus = [3,4,5]

# lines between buses. assume that the graph is connected
const lines = ((1,2), (1,4), (2,3), (2,4), (3,6), (4,5), (5,6))


# if the wind generator is nondispatchable, then treat it like a negative load
if(caseID == "nondispatchable_wind")
	numLoads = 4
	numWindGen = 0
	numGen = 2
	loadToBus = [2,3,4,5]
end



#*========== GENERATOR BOUNDS ===========

# lower and upper power bounds for regular generators (pmin and pmax)
const regGenMinPow = zeros(Float64,numRegGen)
const regGenMaxPow = zeros(Float64,numRegGen)
# set high enough upper bound for the most expensive (second) generator to ensure relatively complete recourse
regGenMaxPow[1] = 220.0
regGenMaxPow[2] = 1000.0

# lower power bounds for wind generators. upper bounds are random
const windGenMinPow = zeros(Float64,numWindGen)



#*========== COST COEFFICIENTS ===========

# regular generator costs including cost of reserves
# only use linear cost coefficients
const regGenCost = [13.5,17.7]

# scale reserves cost by suitable factor
const regGenPosResCost = reservesCostFactor*regGenCost
const regGenNegResCost = 1.0*regGenPosResCost

# wind generator costs including cost of reserves (assume zero marginal costs)
const windGenCost = zeros(Float64,numWindGen)

# set cost of reserves for each wind generator to be a scaled version of minimum regular generator reserves cost
const minRegGenPosResCost = minimum(regGenPosResCost)
const minRegGenNegResCost = minimum(regGenNegResCost)

const windGenPosResCost = minRegGenPosResCost*windReservesCostFactor*ones(Float64,numWindGen)
const windGenNegResCost = minRegGenNegResCost*windReservesCostFactor*ones(Float64,numWindGen)

# allow wind spill (with participation factor 0.1*min_alpha) in the no_wind_reserves case
if(caseID == "no_wind_reserves")
	windGenNegResCost = zeros(Float64,numWindGen)
end



#*========== RESERVES INFORMATION ===========

# which regular generators are ready to provide reserves?
const regGenRes = trues(numRegGen)

# which regular generators are ALWAYS READY to provide reserves?
# make sure that the generator with large enough max generation limit always provides reserves to ensure relatively complete recourse holds
# we usually assume that all generators provide guaranteed reserves with a minimum participation factor min_alpha (that isn't too small) so that numerical issues don't crop up while running the stochastic approximation method
const regGenGuarRes = trues(numRegGen)

# which wind generators are ready to provide reserves?
const windGenRes = trues(numWindGen)
if(caseID == "no_wind_reserves" || caseID == "nondispatchable_wind")
	windGenRes = falses(numWindGen)
end

# which wind generators are ALWAYS READY to provide reserves?
const windGenGuarRes = trues(numWindGen)
if(caseID == "no_wind_reserves" || caseID == "nondispatchable_wind")
	windGenGuarRes = falses(numWindGen)
end

# minimum participation factors for generators ALWAYS READY to provide reserves
# reset this min participation factor to not restrict the space of solutions too much, if necessary
const min_alpha = 0.001
const numGenProvRes = sum(windGenGuarRes) + sum(regGenGuarRes)
if(numGenProvRes > 10)
	min_alpha = 1.0/(100.0*numGenProvRes)
end



#*========== LINE PARAMETERS ===========

# these line parameters are assumed to be in the same order as "lines"
# line susceptances (beta)
const beta = 1./[0.170,0.258,0.037,0.197,0.018,0.037,0.140]

# line flow limits
const lineFlowFactor = 1.0
const flow_max = [200.0,100.0,100.0,100.0,100.0,100.0,100.0]*lineFlowFactor



#*========== LOAD PARAMETERS ===========

# power demands/loads
const demandFactor = 1.0

# average demands
const totalMeanDemand = 280.0*demandFactor
const meanDemand = zeros(Float64,numLoads)
meanDemand[1] = 0.2*totalMeanDemand
meanDemand[2] = 0.4*totalMeanDemand
meanDemand[3] = 0.4*totalMeanDemand



#*========== WIND PARAMETERS ===========

# Wind prediction data
const windFactor = 1.4
const windPredMean = [100.0*windFactor]

# in the nondispatchable case, treat wind generator like a negative load
if(caseID == "nondispatchable_wind")
	meanDemand = [-windPredMean;meanDemand[1:3]]
end



#*========== UNCERTAIN PARAMETERS ===========

# standard deviation of loads/demands
const demandStdevFactor = 0.1
const demandStdev = demandStdevFactor*meanDemand

# randomly generated correlation matrix for the demands
const demandCorrMat = [[1.0 0.5827 0.8654];[0.5827 1.0 0.1761];[0.8654 0.1761 1.0]]

# standard deviation of wind power output
const windPredStdevFactor = 0.1
const windPredStdev = windPredStdevFactor*windPredMean


if(caseID == "nondispatchable_wind")
	# use smaller standard deviation to ensure that relatively complete recourse holds
	demandStdevFactor = 0.05
	demandStdev = demandStdevFactor*meanDemand
	demandCorrMat = [[1.0 0.0 0.0 0.0];[0.0 1.0 0.5827 0.8654];[0.0 0.5827 1.0 0.1761];[0.0 0.8654 0.1761 1.0]]
end



#*========== SCENARIO GENERATION ===========

# store the Cholesky factorization of the demand covariance matrix for numerical efficiency
const demandCovMat = Diagonal(demandStdev)*demandCorrMat*Diagonal(demandStdev)
const demandCovMat_chol = (cholfact(demandCovMat))[:L]


# tailored implementation of multivariate normal distribution for demand scenarios
function generateDemandScenarios(numSamples::Int64)

	omega_d = zeros(Float64,numLoads,numSamples)
	
	omega_d_tmp = rand(Normal(0.0,1.0),numLoads,numSamples)
	for samp = 1:numSamples
		omega_d[:,samp] = demandCovMat_chol*omega_d_tmp[:,samp] + meanDemand
	end
	
	return omega_d
end


# Generate demand and wind scenarios
function generateScenarios(numSamples::Int64)

	omega_d = generateDemandScenarios(numSamples)
	
	omega_w = Array{Float64}(numWindGen,numSamples)
	for i = 1:numWindGen
		if(windPredStdev[i] > 0)
			omega_w[i,:] = windPredMean[i] + rand(Normal(0,windPredStdev[i]), numSamples)		
		else
			omega_w[i,:] = windPredMean[i]*ones(Float64,numSamples)
		end
	end

	return omega_d, omega_w
end