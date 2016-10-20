--[[ This problem is used to valuate a call option. ]]--
volatility = 0.2
interest = 0.00
strike = 20
stockinit = 20
maturity = 1
lower_boundary = math.max(0, math.min(stockinit, strike) - 10)
upper_boundary = math.max(stockinit, strike) + 30

---[[ These parameters define the grid of the problem ]]--
statestep = (upper_boundary - lower_boundary) / 500
statemin = lower_boundary
statemax = upper_boundary
timestep = statestep^2 
	/ (statestep * interest * upper_boundary + (volatility * upper_boundary)^2)

timemax = maturity
policymin = 0.1
policymax = 0.11
policystep = 0.01
tracefile = "result_grid"
optimalvaluefile = "optimal_value"


---[[ The terminal cost of the problem. ]]--
function terminal_cost (state)
	return math.max(state - strike, 0)
end

function boundary_state_max (time, state)
	return math.max((state - statestep) * math.exp(-interest * time) - strike, 0) 
end

function boundary_state_min (time, state)
	return math.max(state * math.exp(-interest * time) - strike, 0) 
end

function running_cost (time, state, control)
	return 0
end

function drift (time, state, control)
	return interest * state
end

function diffusion (time, state, control)
	return volatility * state
end


