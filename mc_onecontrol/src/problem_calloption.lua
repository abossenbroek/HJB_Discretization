--[[ Problem definition ]]--
strike = 20
interest = 0.05
vol = 0.3

---[[ These parameters define the grid of the problem ]]--
statemin = 0
statemax = strike + 50

statestep = math.abs(statemax  - statemin) / 5000

policymin = 0.0001
policymax = 0.0002
policystep = (policymax - policymin) / 1

timestep = 0.001
timemax = 1


--[[ Trace files ]]--
tracefile = "grid"
optimalvaluefile = "optimal_value"
controlfile = "control"

--[[ No running cost for a call option. ]]--
function running_cost (time, control, state)
	return 0
end


---[[ The terminal cost of the problem. ]]--
function terminal_cost (state)
	return math.max(state - strike, 0);
end

function diffusion (time, control, state)
	return vol *  state 
end

function drift (time, control, state)
	return interest * state
end

function boundary_state_max (time, last_state, prev_state)
	--print (time .. " " .. last_state)
	return 2 * last_state - prev_state
end

function boundary_state_min (time, last_state, prev_state)
	return 2 * last_state - prev_state
end


