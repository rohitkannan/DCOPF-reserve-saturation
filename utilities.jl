# UTILITY FUNCTIONS AND DEFINITIONS, INCLUDING FOR CONVERTING USER INPUT TO AMENABLE FORM


#*========== NETWORK TOPOLOGY ===========

# list of regular generators at a particular bus
const regGenAtBus = Array[]

# list of wind generators at a particular bus
const windGenAtBus = Array[]

# list of loads at a particular bus
const loadAtBus = Array[]


# initialize Arrays
for i = 1:numBuses
	push!(regGenAtBus,Int64[])
	push!(windGenAtBus,Int64[])
	push!(loadAtBus,Int64[])
end

# populate regular generators first
for i = 1:numRegGen
	push!(regGenAtBus[regGenToBus[i]],i)
end

# next, populate wind generators
# offset the indices by numRegGen to account for regular generators first
for i = 1:numWindGen
	push!(windGenAtBus[windGenToBus[i]],i+numRegGen)
end

# finally, populate load indices
for j = 1:numLoads
	push!(loadAtBus[loadToBus[j]],j)
end

# store the number of lines
const numLines = size(lines,1)

# store the indices of buses that share an edge/line with current bus
# also store the corresponding line susceptances
const linesAdjBus = Array[]
const betaAdjBus = Array[]

for i = 1:numBuses
	push!(linesAdjBus,Int64[])
	push!(betaAdjBus, Float64[])
end

for k = 1:numLines
	push!(linesAdjBus[lines[k][1]],lines[k][2])
	push!(linesAdjBus[lines[k][2]],lines[k][1])
	
	push!(betaAdjBus[lines[k][1]],beta[k])
	push!(betaAdjBus[lines[k][2]],beta[k])
end


#*========== GENERATOR BOUNDS ===========

# minimum and maximum generator power outputs
const pmin = zeros(Float64,numGen)
const pmax = zeros(Float64,numGen)

# include regular generators first followed by wind generators
pmin[1:numRegGen] = regGenMinPow
pmin[numRegGen+1:end] = windGenMinPow

# ignore the upper bounds for the wind generators since they are random
pmax[1:numRegGen] = regGenMaxPow


# upper bounds for wind generator outputs in the deterministic setting
const pmax_det = 1.0*pmax
for i = 1:numWindGen
	pmax_det[numRegGen+i] = windPredMean[i]
end


# upper bounds for generator outputs to make first-stage decisions
const pmax_FS = 1.0*pmax
for i = 1:numWindGen
	pmax_FS[numRegGen+i] = windPredMean[i] + windBoundFactor*windPredStdev[i]
end


#*========== GENERATOR COST COEFFICIENTS ===========

# generator cost coefficients, including costs of reserves
const genCost = zeros(Float64,numGen)
const posResCost = zeros(Float64,numGen)
const negResCost = zeros(Float64,numGen)

# include regular generators first followed by wind generators
genCost[1:numRegGen] = regGenCost
genCost[numRegGen+1:end] = windGenCost

posResCost[1:numRegGen] = regGenPosResCost
posResCost[numRegGen+1:end] = windGenPosResCost

negResCost[1:numRegGen] = regGenNegResCost
negResCost[numRegGen+1:end] = windGenNegResCost


#*========== GENERATOR RESERVES INFORMATION ===========

# which generators are ready to provide reserves?
const genRes = falses(numGen)

# which generators are ALWAYS READY TO PROVIDE RESERVES?
const genGuarRes = falses(numGen)

# populate the regular generator options followed by the wind generator options
genRes[1:numRegGen] = regGenRes
genRes[numRegGen+1:end] = windGenRes

genGuarRes[1:numRegGen] = regGenGuarRes
genGuarRes[numRegGen+1:end] = windGenGuarRes


#*========== INITIAL GUESSES ===========

# power commitment variables for regular and wind generators
const p0_init = zeros(Float64,numGen)

# positive reserve variables for regular and wind generators
const rp_init = zeros(Float64,numGen)

# negative reserve variables for regular and wind generators
const rm_init = zeros(Float64,numGen)

# participation factor variables for regular and wind generators
const alpha_init = ones(Float64,numGen)/numGen


#*========== OTHER UTILITIES ===========

# update generator upper bounds based on wind levels
function updateGenBounds(omega_w::Array{Float64})

	for i = 1:numWindGen
		if(omega_w[i] < pmin[numRegGen+i])
			println("Wind availability at the ",i,"th generator = ",omega_w[i],"   < corresponding pmin = ",pmin[numRegGen+i])
			throw(ErrorException("Invalid wind generator upper bound!"))
		end
		pmax[numRegGen+i] = omega_w[i]
	end
end

# compute the variance and stdev of sum of demand deviations for the CAP model
const sumDemandDeviationVar = ones(Float64,numLoads)'*demandCovMat*ones(Float64,numLoads)
const sumDemandDeviationStdev = sqrt(sumDemandDeviationVar)


#*========== FACTORIZE AND STORE MATRICES FOR LINEAR SOLVES ===========

# pre-factorize and store matrix corresponding to DC-OPF equations
# the last bus is chosen as the reference bus and substituted out
const rowIndices = []
const colIndices = []
const matValues = []

const sumBeta = zeros(Float64,numBuses-1)

# also construct a matrix of lines with the reference bus ignored
const Line_mat = zeros(Float64,numLines,numBuses-1)


for k = 1:numLines
	fromBus::Int64 = lines[k][1]
	toBus::Int64 = lines[k][2]
	
	if(fromBus < numBuses)
		Line_mat[k,fromBus] = 1.0
		if(toBus < numBuses)
			push!(rowIndices,fromBus)
			push!(colIndices,toBus)
			push!(matValues,-beta[k])
		end
		
		sumBeta[fromBus] += beta[k]
	end
	
	if(toBus < numBuses)
		Line_mat[k,toBus] = -1.0
		if(fromBus < numBuses)
			push!(rowIndices,toBus)
			push!(colIndices,fromBus)
			push!(matValues,-beta[k])
		end
		
		sumBeta[toBus] += beta[k]
	end
end

for i = 1:numBuses-1
	push!(rowIndices,i)
	push!(colIndices,i)
	push!(matValues,sumBeta[i])
end

const B_mat = sparse(rowIndices,colIndices,matValues)
const B_fact = lufact(B_mat)


# construct matrices for generating the rhs of the DCOPF equations
# it appears that storing these in a sparse format is not worth it
# first, construct the matrix corresponding to the generators
const Gen_mat = zeros(Float64,numBuses-1,numGen)

for i = 1:numBuses-1
	for k = 1:size(regGenAtBus[i],1)
		Gen_mat[i,regGenAtBus[i][k]] = 1.0
	end
	
	for k = 1:size(windGenAtBus[i],1)
		Gen_mat[i,windGenAtBus[i][k]] = 1.0
	end
end


# next, construct the matrix corresponding to the loads
const Load_mat = zeros(Float64,numBuses-1,numLoads)

for i = 1:numBuses-1
	for k = 1:size(loadAtBus[i],1)
		Load_mat[i,loadAtBus[i][k]] = 1.0
	end
end

# store the solution of B_mat*x = Gen_mat
const B_over_G = B_fact\Gen_mat

# store the solution of B_mat*x = Load_mat
const B_over_L = B_fact\Load_mat

# store (Line_mat*B_over_G)'
const Line_mat_times_B_over_G_trans = (B_over_G')*(Line_mat')