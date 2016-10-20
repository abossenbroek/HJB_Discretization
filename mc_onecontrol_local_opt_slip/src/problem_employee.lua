--[[ Problem definition ]]--
initial_wealth = 20
risk_aversion = -0.01
expected_return = 0.500
vol = 15
bonus_discount = 0
fixed_wealth = 30

---[[ These parameters define the grid of the problem ]]--
statemin = -10 - fixed_wealth
statemax = initial_wealth + 100

statestep = math.abs(statemax  - statemin) / 1000

policymin = 0.0001
policymax = 2
policystep = (policymax - policymin) / 100

--[[ timestep = statestep^2 
	/ (statestep * expected_return * statemax + (vol * statemax)^2) ]]--
timestep = 0.01

timemax = 1


--[[ Trace files ]]--
tracefile = "grid"
optimalvaluefile = "optimal_value"
controlfile = "control"

--[[ The running cost of the problem is some utility function ]]--
function running_cost (time, control, state)
	bonus = math.max(state + fixed_wealth, 0)
	--[[ The exponential utility function is convex and therefore 
	     models an individual who is risk seeking ]]--
	-- return 1 - math.exp(-risk_aversion * bonus)
	--return math.exp(-risk_aversion * bonus)
	
	--[[ The log utility function is concave and therefore 
	     models an individual who is risk averse ]]--
	return 10 * math.log((bonus + 0.001)) --[[ - control ]]--
	
	--[[ The CRRA (Constant Relative Risk Aversion) utility ]]--
	-- return math.pow(bonus, risk_aversion) / risk_aversion
end


---[[ The terminal cost of the problem. ]]--
function terminal_cost (state)
	return 0
end

function diffusion (time, control, state)
	return vol
end

function drift (time, control, state)
	return control * state * expected_return
end

function boundary_state_max (time, last_state, prev_state)
	--print (time .. " " .. last_state)
	return 2 * last_state - prev_state
end

function boundary_state_min (time, last_state, prev_state)
	return 2 * last_state - prev_state
end


