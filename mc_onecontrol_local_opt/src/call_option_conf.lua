--[[ This problem is used to valuate a call option. ]]--
volatility = 0.2
interest = 0.05
strike = 20
stockinit = 20
maturity = 1
lower_boundary = math.max(1, math.min(stockinit, strike) - 20)
upper_boundary = math.max(stockinit, strike) + 80

---[[ These parameters define the grid of the problem ]]--
statestep = (upper_boundary - lower_boundary) / 400
statemin = lower_boundary + statestep
statemax = upper_boundary
 timestep = statestep^2 
 	/ (statestep * interest * upper_boundary + (volatility * upper_boundary)^2)   
--timestep = 0.1
-- timestep = 0.0225


timemax = maturity
policymin = 0.1
policymax = 0.11
policystep = 0.01
tracefile = "result_grid"
optimalvaluefile = "call_option_optimal_value"
controlfile = "opt_control"
