---[[ These parameters define the grid of the problem ]]--
statestep = 0.1
statemin = 1
statemax = 20

timestep = 0.1
timemax = 1


policymin = 0.1
policymax = 1
policystep = 0.02


--[[ Trace files ]]--
tracefile = "grid"
optimalvaluefile = "optimal_value"
controlfile = "control"


---[[ The terminal cost of the problem. ]]--
function terminal_cost (state)
	return 1;
end

function boundary_state_max (time, last_state, prev_state)

end

function boundary_state_min (time, last_state, prev_state)
end

function running_cost (time, control, state)
end

function drift (time, state, control)
end

function diffusion (time, control, state)
end


