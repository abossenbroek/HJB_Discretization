--[[ Problem definition ]]--
strike = 20
vol = 0.2
interest = 0.20

---[[ These parameters define the grid of the problem ]]--
statemin = 0
statemax = strike + 50

statestep = math.abs(statemax  - statemin) / 1000

policymin = 0.0001
policymax = 0.0002
policystep = (policymax - policymin) / 1

timestep = 0.01
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
	return 50
end

function boundary_state_min (time, last_state, prev_state)
	return 0
end


