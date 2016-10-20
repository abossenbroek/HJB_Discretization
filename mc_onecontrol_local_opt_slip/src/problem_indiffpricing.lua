--[[ Problem definition ]]--
mu = 0.2
r = 0
sigma = 0.6
initial_wealth = 00
gamma = 0.5

---[[ These parameters define the grid of the problem ]]--
statemin = 0
statemax = initial_wealth + 50

statestep = math.abs(statemax - statemin) / 1000

policymin = 0.1
policymax = 20
policystep = (policymax - policymin) / 250

--timestep = statestep^2 
--	/ (statestep * mu * policymax * statemax + (sigma * policymax * statemax)^2) 
timestep = 0.005
timemax = 1


--[[ Trace files ]]--
tracefile = "grid"
optimalvaluefile = "optimal_value"
controlfile = "control"

--[[ The running cost of the problem is some utility function ]]--
function running_cost (time, control, state)
	--[[ The exponential utility function is convex and therefore 
	     models an individual who is risk seeking ]]--
	return -math.exp(-gamma * state)
	
	--[[ The log utility function is concave and therefore 
	     models an individual who is risk averse ]]--
	-- return math.log(gamma * state) 
	
	--[[ The CRRA (Constant Relative Risk Aversion) utility ]]--
	-- return math.pow(state, gamma) / gamma
end


---[[ The terminal cost of the problem. ]]--
function terminal_cost (state)
	return -math.exp(-gamma * state)
end

function diffusion (time, control, state)
	return sigma * control
end

function drift (time, control, state)
	return (r + control * (mu - r))
end

function boundary_state_max (time, last_state, prev_state)
	--print (time .. " " .. last_state)
	return 2 * last_state - prev_state
end

function boundary_state_min (time, last_state, prev_state)
	return 2 * last_state - prev_state
end


