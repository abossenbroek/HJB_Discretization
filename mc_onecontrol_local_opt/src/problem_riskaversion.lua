--[[ Problem definition ]]--
initial_wealth = 20
risk_aversion = -0.2
expected_return = 0.100
vol = 0.2

---[[ These parameters define the grid of the problem ]]--
statemin = 0.1
statemax = initial_wealth + 40

statestep = math.abs(statemax  - statemin) / 200

policymin = 0.0001
policymax = 1
policystep = (policymax - policymin) / 200

timestep = statestep^2 
	/ (statestep * expected_return * statemax + (vol * statemax)^2)
--timestep = 0.01

timemax = 1


--[[ Trace files ]]--
tracefile = "grid"
optimalvaluefile = "optimal_value"
controlfile = "control"

--[[ The running cost of the problem is some utility function ]]--
function running_cost (time, control, state)
	--[[ The exponential utility function is convex and therefore 
	     models an individual who is risk seeking ]]--
	return math.exp(-risk_aversion * state)
	--[[ The log utility function is concave and therefore 
	     models an individual who is risk averse ]]--
	--return math.log(risk_aversion * state)
end


---[[ The terminal cost of the problem. ]]--
function terminal_cost (state)
	return 0
end

function diffusion (time, control, state)
	return (1 - control) * state * vol
end

function drift (time, control, state)
	--print ("drift for state " .. state .. " and control " .. control .. " is " .. control * state * expected_return)
	return control * state * expected_return
end

function boundary_state_max (time, last_state, prev_state)
	--print (time .. " " .. last_state)
	return 2 * last_state - prev_state
end

function boundary_state_min (time, last_state, prev_state)
	return 2 * last_state - prev_state
end


