# 118 Bus Model based on the IEEE 118 bus case data in http://motor.ece.iit.edu/data/JEAS_IEEE118.doc
# Wind data and modifications partly based on Roald et al's paper "Optimal power flow with wind power control and limited expected risk of overloads"

# OPTIONS for caseID: 
# wind_penetration_level_1: 25% wind penetration
# wind_penetration_level_2: 50% wind penetration
# wind_penetration_level_3: 75% wind penetration
# wind_penetration_level_4: 100% wind penetration
# wind_penetration_level_5: 125% wind penetration

const caseID = "wind_penetration_level_1"



# include basic model parameters
include("params.jl")



#*========== BASIC TOPOLOGY ===========

# number of buses, generators (regular and wind), loads
const numBuses = 118
const numGen = 79
const numRegGen = 54
const numWindGen = 25
const numLoads = 91


# map regular generators to buses
# from column 2 of Table 1 in the reference
const regGenToBus = [4,6,8,10,12,15,18,19,24,25,26,27,31,32,34,36,40,42,46,49,54,55,56,59,61,62,65,66,69,70,72,73,74,76,77,80,82,85,87,89,90,91,92,99,100,103,104,105,107,110,111,112,113,116]

# map wind generators to buses
# from Roald et al's spreadsheet
const windGenToBus = [3,6,9,18,25,29,33,35,41,43,46,56,60,61,70,78,83,85,92,93,96,107,113,117,118]

# map loads to buses
# from column 1 of Table 6 in the reference
const loadToBus = [1,2,3,4,6,7,11,12,13,14,15,16,17,18,19,20,21,22,23,27,28,29,31,32,33,34,35,36,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,62,66,67,70,74,75,76,77,78,79,80,82,83,84,85,86,88,90,92,93,94,95,96,97,98,100,101,102,103,104,105,106,107,108,109,110,112,114,115,117,118]

# lines between buses
# from columns 2 and 3 of Table 3 in the reference
const lines = ((1,2), (1,3), (4,5), (3,5), (5,6), (6,7), (8,9), (8,5), (9,10), (4,11), (5,11), (11,12), (2,12), (3,12), (7,12), (11,13), (12,14), (13,15), (14,15), (12,16), (15,17), (16,17), (17,18), (18,19), (19,20), (15,19), (20,21), (21,22), (22,23), (23,24), (23,25), (26,25), (25,27), (27,28), (28,29), (30,17), (8,30), (26,30), (17,31), (29,31), (23,32), (31,32), (27,32), (15,33), (19,34), (35,36), (35,37), (33,37), (34,36), (34,37), (38,37), (37,39), (37,40), (30,38), (39,40), (40,41), (40,42), (41,42), (43,44), (34,43), (44,45), (45,46), (46,47), (46,48), (47,49), (42,49), (42,49), (45,49), (48,49), (49,50), (49,51), (51,52), (52,53), (53,54), (49,54), (49,54), (54,55), (54,56), (55,56), (56,57), (50,57), (56,58), (51,58), (54,59), (56,59), (56,59), (55,59), (59,60), (59,61), (60,61), (60,62), (61,62), (63,59), (63,64), (64,61), (38,65), (64,65), (49,66), (49,66), (62,66), (62,67), (65,66), (66,67), (65,68), (47,69), (49,69), (68,69), (69,70), (24,70), (70,71), (24,72), (71,72), (71,73), (70,74), (70,75), (69,75), (74,75), (76,77), (69,77), (75,77), (77,78), (78,79), (77,80), (77,80), (79,80), (68,81), (81,80), (77,82), (82,83), (83,84), (83,85), (84,85), (85,86), (86,87), (85,88), (85,89), (88,89), (89,90), (89,90), (90,91), (89,92), (89,92), (91,92), (92,93), (92,94), (93,94), (94,95), (80,96), (82,96), (94,96), (80,97), (80,98), (80,99), (92,100), (94,100), (95,96), (96,97), (98,100), (99,100), (100,101), (92,102), (101,102), (100,103), (100,104), (103,104), (103,105), (100,106), (104,105), (105,106), (105,107), (105,108), (106,107), (108,109), (103,110), (109,110), (110,111), (110,112), (17,113), (32,113), (32,114), (27,115), (114,115), (68,116), (12,117), (75,118), (76,118))



#*========== GENERATOR BOUNDS ===========

# lower and upper power bounds for regular generators (pmin and pmax)
const regGenMinPow = zeros(Float64,numRegGen)
# from column 6 of Table 1 in the reference
# change upper limit of generator#49 to be large enough to ensure relatively complete recourse holds
const regGenMaxPow = [30.0,30.0,30.0,300.0,300.0,30.0,100.0,30.0,30.0,300.0,350.0,30.0,30.0,100.0,30.0,100.0,30.0,30.0,100.0,250.0,250.0,100.0,100.0,200.0,200.0,100.0,420.0,420.0,300.0,80.0,30.0,30.0,20.0,100.0,100.0,300.0,100.0,30.0,300.0,200.0,20.0,50.0,300.0,300.0,300.0,20.0,100.0,100.0,10000.0,50.0,100.0,100.0,100.0,50.0]

# lower power bounds for wind generators
const windGenMinPow = zeros(Float64,numWindGen)



#*========== COST COEFFICIENTS ===========

# regular generator costs including cost of reserves
# only use linear cost coefficients
# from column 4 of Table 1 in the reference
const regGenCost = [26.2438,26.2438,26.2438,12.8875,12.8875,26.2438,17.8200,26.2438,26.2438,12.8875,10.7600,26.2438,26.2438,17.8200,26.2438,17.8200,26.2438,26.2438,17.8200,12.3299,12.3299,17.8200,17.8200,13.2900,13.2900,17.8200,8.3391,8.3391,12.8875,15.4708,26.2438,26.2438,37.6968,17.8200,17.8200,12.8875,17.8200,26.2438,10.7600,12.8875,37.6968,22.9423,12.8875,12.8875,12.8875,37.6968,17.8200,17.8200,37.6968,22.9423,17.8200,17.8200,17.8200,22.9423]

# scale reserves cost by suitable factor
const regGenPosResCost = reservesCostFactor*regGenCost
const regGenNegResCost = 1.0*regGenPosResCost

# wind generator costs including cost of reserves (assume zero marginal costs)
const windGenCost = zeros(Float64,numWindGen)

# set cost of reserves for each wind generator to be the minimum reserves cost of all regular generators
const minRegGenPosResCost = minimum(regGenPosResCost)
const minRegGenNegResCost = minimum(regGenNegResCost)

const windGenPosResCost = minRegGenPosResCost*windReservesCostFactor*ones(Float64,numWindGen)
const windGenNegResCost = minRegGenNegResCost*windReservesCostFactor*ones(Float64,numWindGen)



#*========== RESERVES INFORMATION ===========

# which regular generators are ready to provide reserves?
const regGenRes = trues(numRegGen)

# which regular generators are ALWAYS READY to provide reserves?
# make sure that the generator with large enough max generation limit always provides reserves to ensure relatively complete recourse holds
# we assume that all generators provide guaranteed reserves with a minimum participation factor min_alpha so that numerical issues don't crop up while running the stochastic approximation method
const regGenGuarRes = trues(numRegGen)

# which wind generators are ready to provide reserves?
const windGenRes = trues(numWindGen)

# which wind generators are ALWAYS READY to provide reserves?
const windGenGuarRes = trues(numWindGen)

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
# from column 6 of Table 3 in the reference
const beta = 1./[0.09990,0.04240,0.00798,0.10800,0.05400,0.02080,0.03050,0.02670,0.03220,0.06880,0.06820,0.01960,0.06160,0.16000,0.03400,0.07310,0.07070,0.24440,0.19500,0.08340,0.04370,0.18010,0.05050,0.04930,0.11700,0.03940,0.08490,0.09700,0.15900,0.04920,0.08000,0.03820,0.16300,0.08550,0.09430,0.03880,0.05040,0.08600,0.15630,0.03310,0.11530,0.09850,0.07550,0.12440,0.24700,0.01020,0.04970,0.14200,0.02680,0.00940,0.03750,0.10600,0.16800,0.05400,0.06050,0.04870,0.18300,0.13500,0.24540,0.16810,0.09010,0.13560,0.12700,0.18900,0.06250,0.32300,0.32300,0.18600,0.05050,0.07520,0.13700,0.05880,0.16350,0.12200,0.28900,0.29100,0.07070,0.00955,0.01510,0.09660,0.13400,0.09660,0.07190,0.22930,0.25100,0.23900,0.21580,0.14500,0.15000,0.01350,0.05610,0.03760,0.03860,0.02000,0.02680,0.09860,0.03020,0.09190,0.09190,0.21800,0.11700,0.03700,0.10150,0.01600,0.27780,0.32400,0.03700,0.12700,0.41150,0.03550,0.19600,0.18000,0.04540,0.13230,0.14100,0.12200,0.04060,0.14800,0.10100,0.19990,0.01240,0.02440,0.04850,0.10500,0.07040,0.02020,0.03700,0.08530,0.03665,0.13200,0.14800,0.06410,0.12300,0.20740,0.10200,0.17300,0.07120,0.18800,0.09970,0.08360,0.05050,0.15810,0.12720,0.08480,0.15800,0.07320,0.04340,0.18200,0.05300,0.08690,0.09340,0.10800,0.20600,0.29500,0.05800,0.05470,0.08850,0.17900,0.08130,0.12620,0.05590,0.11200,0.05250,0.20400,0.15840,0.16250,0.22900,0.03780,0.05470,0.18300,0.07030,0.18300,0.02880,0.18130,0.07620,0.07550,0.06400,0.03010,0.20300,0.06120,0.07410,0.01040,0.00405,0.14000,0.04810,0.05440]

# line flow limits
const lineFlowFactor = 0.75
# from column 8 of Table 3 in the reference
const flow_max = [175.0,175.0,500.0,175.0,175.0,175.0,500.0,500.0,500.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,500.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,500.0,500.0,500.0,175.0,175.0,500.0,175.0,500.0,175.0,175.0,140.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,500.0,500.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,500.0,175.0,175.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,175.0,175.0,500.0,175.0,500.0,175.0,175.0,500.0,500.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,500.0,175.0,175.0,175.0,175.0,175.0,175.0,500.0,500.0,175.0,500.0,500.0,200.0,200.0,175.0,175.0,175.0,500.0,500.0,175.0,175.0,500.0,500.0,500.0,175.0,500.0,500.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,200.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,500.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,500.0,175.0,175.0,175.0,500.0,175.0,175.0,175.0]*lineFlowFactor



#*========== LOAD PARAMETERS ===========

# power demands/loads
const demandFactor = 1.5

# from column 2 of Table 6 in the reference
const meanDemand = [54.14,21.23,41.40,31.85,55.20,20.17,74.31,49.89,36.09,14.86,95.54,26.54,11.68,63.69,47.77,19.11,14.86,10.62,7.43,65.82,18.05,25.48,45.65,62.63,24.42,62.63,35.03,32.91,27.00,20.00,37.00,37.00,18.00,16.00,53.00,28.00,34.00,20.00,87.00,17.00,17.00,18.00,23.00,113.00,63.00,84.00,12.00,12.00,277.00,78.00,77.00,39.00,28.00,66.00,68.00,47.00,68.00,61.00,71.00,39.00,130.00,54.00,20.00,11.00,24.00,21.00,48.00,78.00,65.00,12.00,30.00,42.00,38.00,15.00,34.00,37.00,22.00,5.00,23.00,38.00,31.00,43.00,28.00,2.00,8.00,39.00,25.00,8.49,23.35,21.23,33.00]*demandFactor



#*========== WIND PARAMETERS ===========

# Wind prediction data
const windFactor = 1.0

# from Roald et al's spreadsheet
const windPredMean = [109.0,70.0,200.0,147.0,102.0,148.0,70.0,60.0,96.0,80.0,105.0,113.0,80.0,84.0,59.0,146.0,69.0,250.0,118.0,89.0,127.0,76.0,72.0,129.0,129.0]

const sumMeanWind = sum(windPredMean)
const sumMeanDemand = sum(meanDemand)

# set average wind output depending on the wind penetration level
if(caseID == "wind_penetration_level_1")
	windFactor = 0.25*sumMeanDemand/sumMeanWind
elseif(caseID == "wind_penetration_level_2")
	windFactor = 0.5*sumMeanDemand/sumMeanWind
elseif(caseID == "wind_penetration_level_3")
	windFactor = 0.75*sumMeanDemand/sumMeanWind
elseif(caseID == "wind_penetration_level_4")
	windFactor = 1.0*sumMeanDemand/sumMeanWind
elseif(caseID == "wind_penetration_level_5")
	windFactor = 1.25*sumMeanDemand/sumMeanWind
end

windPredMean *= windFactor



#*========== UNCERTAIN PARAMETERS ===========

# standard deviation of loads/demands. keep this small enough to be realistic
const demandStdevFactor = 0.05
const demandStdev = demandStdevFactor*meanDemand

# randomly generated correlation matrix for the demands
include("118bus_demandCorrMat.jl")

# standard deviation of wind power output
const windPredStdevFactor = 0.1
const windPredStdev = windPredStdevFactor*windPredMean

# randomly generated correlation matrix for the wind power outputs
include("118bus_windCorrMat.jl")



#*========== SCENARIO GENERATION ===========

# store the Cholesky factorization of the demand covariance matrix for numerical efficiency
const demandCovMat = Diagonal(demandStdev)*demandCorrMat*Diagonal(demandStdev)
const demandCovMat_chol = (cholfact(Hermitian(demandCovMat)))[:L]

# store the Cholesky factorization of the wind covariance matrix for numerical efficiency
const windCovMat = Diagonal(windPredStdev)*windCorrMat*Diagonal(windPredStdev)
const windCovMat_chol = (cholfact(Hermitian(windCovMat)))[:L]


# tailored implementation of multivariate normal distribution for demand scenarios
function generateDemandScenarios(numSamples::Int64)

	omega_d = zeros(Float64,numLoads,numSamples)
	
	omega_d_tmp = rand(Normal(0.0,1.0),numLoads,numSamples)
	for samp = 1:numSamples
		omega_d[:,samp] = demandCovMat_chol*omega_d_tmp[:,samp] + meanDemand
	end
	
	return omega_d
end


# tailored implementation of multivariate normal distribution for wind scenarios
function generateWindScenarios(numSamples::Int64)

	omega_w = zeros(Float64,numWindGen,numSamples)
	
	omega_w_tmp = rand(Normal(0.0,1.0),numWindGen,numSamples)
	for samp = 1:numSamples
		omega_w[:,samp] = windCovMat_chol*omega_w_tmp[:,samp] + windPredMean
	end
	
	return omega_w
end



# Generate demand and wind scenarios
function generateScenarios(numSamples::Int64)

	omega_d = generateDemandScenarios(numSamples)
	omega_w = generateWindScenarios(numSamples)

	return omega_d, omega_w
end