### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ f59fce3a-3273-11ec-3fe8-9fed7d8c159e
begin
	using Random
	using Plots
	using Statistics, StatsBase
	using Plots.PlotMeasures
	using PlutoUI
	using UnPack
	using ThreadsX
	import Base.show
	using GLM
	using DataFrames
	using Optim
end; md"### Import packages"

# ╔═╡ 74fb6a63-48a4-4f52-bf49-dfcc70785f66
md"""
# Dynamical models of the risky choice behavior and perturbations

We found that bilateral muscimol infusion into the FOF dramatically shifted choices towards the surebet. Bilateral optogenetic silencing gave a consistent result.

In contrast, unilateral silencing with muscimol had a very small **non-lateralized** effect. 

In our Biorxiv preprint, we included a single dynamical model as an existence proof - to show how bilateral silencing might give the observed results.  To be more systematic we have added several additional models and also added a decision process to our models, so that we can make psychometric curves of the control and perturbed models.

"""

# ╔═╡ 0a41f9e0-36d8-406e-8fca-f6f277b03d1c
md"""
### FOF as 1/6 of a distributed network. 

Our first model was inspired by the results from Kopec (2015), that the FOF is 1/6 of a distributed network. However, this model differs from that one in that we are assuming that the unilateral effects are negligble (as in our data) and that when we silence the FOF we are simulating a bilateral inhibition. This simulation was originally designed to explain our muscimol results and predicted that neurons in FOF would encode the lottery value.

By movind the `seed` slider, you can test the simulation with different randon weights.
"""

# ╔═╡ 0d028ece-4112-4518-b320-f896704c97ab
md"""
seed = $(@bind Σ Slider(1:1000, show_value = true, default=205))

$(@bind print_it Button("Save Figure"))
"""

# ╔═╡ 02ac31e6-e9d1-4946-bbd1-9ed802310d63
md"""
### The FOF as a two-node network with left and right driving left/right choices.

This is basically reimplementation of the model in Kopec 2015.

The prediction here is that bilateral inactivation would just make choices noisier, it would also predict very different outcomes of unilateral inactivation and of recordings
"""

# ╔═╡ eec8c196-9808-44db-a65b-dc3c52b8c04e
 

# ╔═╡ d38ef57a-2519-4f95-b5a2-8d890b29915e
md"""
## Code Appendix
"""

# ╔═╡ c28c212e-3f37-4ad5-ae4b-9377784a22f0
begin
	"""
	Options

	A struct to hold options for the simulation. Could use a dictionary instead, but i like the . notation
	"""
	mutable struct Options
		n_units::Int16
		w_e::Float32
		w_e_σ::Float32
		tau::Float32
		dt::Float32
	    starttime::Float32
		endtime::Float32
		raw_σ::Float32
		scale::Float32
		σ::Float32
		timex::AbstractRange
	end
	Options(n_units, w_e,w_e_σ, tau, dt, starttime, endtime, raw_σ, scale) = 
	Options(n_units, w_e, w_e_σ,tau, dt, starttime, endtime, raw_σ, scale, (raw_σ/dt)/n_units,starttime:dt:endtime)
	# Here is where we set the default options.
	# OLD DEFAULTS FOR TAU INSIDE.
	# Options() = Options(6,5,1,0.15,0.001, -5, 1.5, 3)
	Options() = Options(6,5*0.15,.1,0.15,0.001, -5, 1.5, 2*0.15, 1.0)
	OptionsLR() = Options(2, 5*0.15,.01,0.15,0.001, -5, 1.5, 0.15, 1.0)
	Base.show(io::IO,::MIME"text/html",x::Options) = begin
		for f in fieldnames(typeof(x))
	 		write(io, "<b>$f</b>: $(getfield(x,f))<br>")
		end
	end

end; md"### Options"

# ╔═╡ 40f6aa19-583b-4e88-93e7-fabf16fd082f
"""
ci(x)

returns 95% CI in a format that is easily used. 
```
plot(x,mean(y),yerr=ci(y))
```
"""
ci(x::Matrix, boots=200) = begin
		mu = mean.(eachcol(x))
		#q = std(eachcol(x))
		#(mu - q, mu+q)
		mu_ci = zeros(eltype(x), boots, length(mu))
		for bx ∈ 1:boots
			mu_ci[bx,:] = mean(hcat(
				[sample(z, length(z)) for z in eachcol(x)]...)
				, dims=1)
		end
		q = quantile.(eachcol(mu_ci), Ref([0.025, 0.975]))
		(abs.(mu .- first.(q)), abs.(last.(q) .- mu))
	
	end

# ╔═╡ f4ba91f5-aba0-46a3-a8b6-b2da0fb4bbea
md"""### Code for original model"""

# ╔═╡ 917f1cd3-6545-42de-a2ea-e48b9d6f9665
"""
fit_sim(T)

T is the 2nd output of `do_sim`

Estimate the ρ in U=V^ρ
"""
fit_sim(T) = begin
 x =  vcat(repeat(T[1],inner=(20,1))...)
 cntrl =  vcat(T[2]...)
 opto  =  vcat(T[3]...)
	mse(ρ, X) = begin
		y_hat = X[1] .^ ρ[1] .+ ρ[2]
		return sum( (y_hat .- X[2]).^2 )

	end

	x0 = [1.0,1.0]
	f(z) = mse(z, (x,cntrl))
	o1 = optimize(f, x0, BFGS())
	g(z) = mse(z, (x,opto))
	o2 = optimize(g, x0, BFGS())
	return o1, o2
	
end

# ╔═╡ 8d03feb6-8d94-4765-914a-266b7b04034e
logit(params,x) = begin
	x_0, b = params
	return 1 ./ (1 .+ exp.(-1 .* b .* (x .- x_0)))
end

# ╔═╡ 903f6ce1-66e8-470b-815f-2a6a69581d9a
"""
fitlogit(p, x, y)

for some parameters `p` return the MSE of the choices generates by `p` compared to `y` given trials `x`.
"""
fitlogit(params,x,y) = begin
	x_0, b = params
	p = logit(params,x)
	ll(a,b) = (a-b)^2
	return sum(ll.(y,p))
end

# ╔═╡ 3610a639-8a8c-4af6-b161-0c5badfc3e8e
"""
psycho_sim(T)

generate a psychometric curve based on simulated data returned from do_sim (or something similar) 
"""
psycho_sim(T; β=0.2, x₀=NaN, as_choice=false) = begin
	lotteries, control, muscimol = T
	if isnan(x₀)
		x₀ = quantile(control[:],0.5)
		# Make it so they choose lottery 65% of the time in the cntrl simulation
	end
	softmax(x) = 1 / (1 + exp(-β * (x - x₀)))
	if as_choice
		choose(x) = (rand() < 1 / (1 + exp(-β * (x - x₀))))
	else
		choose = softmax
	end
	X = reshape(repeat(lotteries, inner=size(control,1)), 20,6)
	Y_control = choose.(control)
	Y_musc = choose.(muscimol)
	#lm = glm(@formula(chose_lottery ~ value), DataFrame(chose_lottery=Y_control[:], value=X[:]),Bernoulli)

	f1(p) = fitlogit(p, X[:] .+ 0.0, Y_control[:])
	o1 = optimize(f1, [51.0 ,0.1])
	p1 = o1.minimizer

	f2(p) = fitlogit(p, X[:] .+ 0.0, Y_musc[:])
	o2 = optimize(f2, [51.0 ,0.1])
	p2 = o2.minimizer
	
	
	return (X, Y_control, Y_musc, p1, p2)
	
	
	
end
	
		
		

# ╔═╡ f0e68369-c342-43e1-80d3-71ef27eede36
md"""#### defaults parameters for model1

$(Options())

"""

# ╔═╡ 3aee647e-6150-448f-af48-6ddcadd4f688
"""
simMint(opto; opts = Options(), W=[], input=[],value=0)

`opto` is the fraction of units to silence. 

Simulates the dynamical system using Euler's method with defaults set by `opts`. Outputs the timecourse of the simulation for `trace_plot`
"""
simMint(opto; opts = Options(), W=[], input=[],value=0) = begin

	
	@unpack n_units, w_e, tau, dt, timex, σ, w_e_σ, scale = opts

    U = fill(0.0, (length(timex), n_units))
	#W = ones(eltype(U), (n_units,n_units))*w_e/n_units # try random?
	if isempty(W)
		W = randn(n_units,n_units)*w_e_σ .+ w_e/n_units
	end
	if isempty(input)
		input = copy(U)
		input[0.1 .< timex .< 1.5,floor(Int16,end/2):end] .= value^scale * tau
	end	
	opto_units = [x/ n_units < opto for x in 1:n_units]    
	
	for t = 2:length(timex)

		# @show size.([U[t,:], U[t-1,:], input_on, input_off] )
		
		U[t,:] = U[t - 1,:] + 
			(dt / tau) * (-(U[t - 1,:]) +
			      W*U[t - 1,:] + 
				 input[t,:] + randn(n_units) * σ)

		U[t, U[t,:].<0] .= 0.0
		U[t, U[t,:].>100] .= 100.0
		U[t,opto_units] .= 0.0

    end
    return (timex, mean.(eachrow(U)),U, W) 
end

# ╔═╡ 10306590-a806-4708-a823-07d11d5c1c2d
"""
simM(value, opto; opts = Options(), W=[])

A wrapper for `simMint` (just to save me some work because I restructured the inputs)
"""
simM(value, opto; opts = Options(), W=[]) = begin

	simMint(opto; value, opts, W)
end

# ╔═╡ 2263c34f-1481-40f0-9893-32334ce65009
begin
	my_seed = Σ
	Random.seed!(my_seed)
	t,μ,U,W = simM(12,0)
end; md"""#### Get a good weight matrix"""

# ╔═╡ 5789c252-9ffa-479c-8bd1-bad310ee9005
"""
simMfast(value::Number, opto::Number; opts = Options(), W = [])

`value` is the value of the lottery. `opto` is the fraction of the network to silence.
"""
simMfast(value::Number, opto::Number; opts = Options(), W = []) = begin

	simMfast(opto; value, opts, W)
end; 

# ╔═╡ 752ea43e-5f42-410c-acfe-2dcc795b4c89
"""
simMfast(opto::Number; opts = Options(), W=[], input=[],value=0.0)

a faster network simluation. updates the network nodes inplace to save on memory.
"""
simMfast(opto::Number; opts = Options(), W=[], input=[],value=0.0) = begin

	
	@unpack n_units, w_e, tau, dt, timex, σ, w_e_σ, scale = opts

    U = fill(0.0, n_units)
	#W = ones(eltype(U), (n_units,n_units))*w_e/n_units # try random?
	if isempty(W)
		W = randn(n_units,n_units)*w_e_σ .+ w_e/n_units
	end
	if isempty(input)
		input = fill(0.0, (length(timex), n_units))
		input[0.1 .< timex .< 1.5,floor(Int16,end/2):end] .= value^scale * tau 
	end	
	opto_units = [x/ n_units < opto for x in 1:n_units]    
	
	for t = 2:length(timex)

		# @show size.([U[t,:], U[t-1,:], input_on, input_off] )
		
		U += (dt/tau) * (-U +
			      W*U + 
				 input[t,:] + randn(n_units) * σ)

		U[U.<0] .= 0.0
		U[U.>100] .= 100.0
		U[opto_units] .= 0.0

    end
    return (timex, mean(U),U, W) 
end

# ╔═╡ 030cea43-fecf-417a-9e1d-7a6970569234
"""
`do_sim`

Wrapper for generating the FOF as 1/6 simulation
"""
do_sim(;boots=20, opts=Options(), W=[]) = begin
	values = [12*x for x in [0,1,2,4,8,16]]
	cntrl = fill(0.0, boots, length(values))
	opto = similar(cntrl)
	for (vi,value) = enumerate(values)
			cntrl[:, vi] = ThreadsX.map(x->simMfast(value, 0; opts,W)[2], 1:boots)
			opto[:, vi] = ThreadsX.map(x->simMfast(value, 0.2; opts,W)[2], 1:boots)

	end
	return (values, cntrl, opto)
end;md"#### do_sim"

# ╔═╡ a4a67fe2-3dcd-4d0f-a514-d15ab736e01e
begin
	T = do_sim(boots=20;W); md"#### Run Simulation"
	sb_val = quantile(T[2][:],0.65)
	o1,o2 = fit_sim(T)
end; md"#### Run simulation"

# ╔═╡ aaebd031-a8de-40d8-a2da-8c9408d89e0d
"""
make_input

"""
make_input(value, step, opts; 
	starttime=0.1, endtime=1.5, steptime=0.7, input_frac=0.5) = begin
t = opts.timex
inp_units = [i/opts.n_units > input_frac for i ∈ 1:opts.n_units]
I = fill(0.0, (length(t), opts.n_units))
I[starttime.<=t.<steptime,inp_units] .= value + step
I[steptime.<=t.<endtime,inp_units] .= value - step
I		
end#; md"#### make input"

# ╔═╡ 923cacc6-f84d-44c3-829b-9b82d121a80c
"""
do_sim_leaky()

I'm not totally sure why i wrote this. It tests how early vs. late input influences the network. 
"""
do_sim_leaky(;boots=20, opts=Options(), W=[]) = begin
	
	timex = opts.timex
	inputH(value) = make_input(value, 30, opts)
	inputL(value) = make_input(value, -30, opts)
	values = [12*x for x in [0,1,2,4,8,16]]
	cntrl_high = fill(0.0, boots, length(values))
	cntrl_low = similar(cntrl_high)
	opto_high = similar(cntrl_high)
	opto_low = similar(cntrl_high)

	for (vi,value) = enumerate(values)
		cntrl_high[:, vi] = 
			ThreadsX.map(x->simMfast(0; opts,W, input=inputH(value))[2], 1:boots);
		cntrl_low[:, vi] = 
			ThreadsX.map(x->simMfast(0; opts,W,input=inputL(value))[2], 1:boots);
		opto_high[:, vi] = 
			ThreadsX.map(x->simMfast(0.2; opts,W,input=inputH(value))[2], 1:boots);
		opto_low[:, vi] = 
			ThreadsX.map(x->simMfast(0.2; opts,W,input=inputL(value))[2], 1:boots);
	end
	return (values, cntrl_low, cntrl_high, opto_low, opto_high)
end;md"#### do sim leaky"

# ╔═╡ 65364d78-2b7e-4fc6-80e8-d4a1af740676
md"""
### Code for 2-node L/R model
"""

# ╔═╡ 8eb705a5-f23e-4c60-bfff-5fd2752b363e
"""
simM2trace(opto; opts = OptionsLR(), W=[1 -1; -1 1], input=[],value=0)

Simulated the 2-node L/R model

`opto` is a tuple (index, gain) so U[index] *= gain

If `input` is not specified then `value` is the input for the first node and the second node gets constant input of 36. This means that the first node should be considered the lottery node and the 2nd node the surebet node.

Should the input be on all the time or just while the sound is on? for now **yes**. interesting question. this could potentially have a big impact on the results.

Simulates the dynamical system using Euler's method with defaults set by `opts`. Outputs the timecourse of the simulation for `trace_plot`
"""
simM2trace(opto; opts = OptionsLR(), W=[0 -4; -4 0], input=[],value=0) = begin

	@unpack n_units, w_e, tau, dt, timex, σ, w_e_σ, scale = opts
	n_units = 2 # ALWAYS 2 units
    U = fill(0.0, (length(timex), n_units))
	E = [5.0; 5.0]
	#W = ones(eltype(U), (n_units,n_units))*w_e/n_units # try random?
	if isempty(W)
		W = randn(n_units,n_units)*w_e_σ .+ w_e/n_units
	end
	if isempty(input)
		input = copy(U)
		input[0.1 .< timex .< 1.5, 1] .= value^scale * tau
		input[0.1 .< timex .< 1.5, 2] .= 36 * tau
		input[0.1 .< timex .< 1.5, opto[1]] *= opto[2]

	end	
	opto_units = opto[1]    
	
	for t = 2:length(timex)

		# @show size.([U[t,:], U[t-1,:], input_on, input_off] )
		
		U[t,:] = U[t - 1,:] + 
			(dt / tau) * (-(U[t - 1,:]) +
			      W*U[t - 1,:] + E + 
				 input[t,:] + randn(n_units) * σ/3)

		U[t, U[t,:].<0] .= 0.0
		U[t, U[t,:].>80] .= 80.0

    end
    return (timex, mean.(eachrow(U)),U, W) 
end

# ╔═╡ 1d488acb-f9d3-4de8-b193-2016159e8ea2
begin
	WLR = [0 -4; -4 0]
	olr = OptionsLR()
	olr.:raw_σ = 0.1
end; md"#### Play with simulation options"

# ╔═╡ c28e9c6b-5ffa-40f5-be7b-642bca8d6c65
plot = Plots.plot

# ╔═╡ 90b5837b-c763-47ac-91a8-1354c3551e2d


# ╔═╡ 05e638d3-6b71-44f4-a87d-d920bee04b1b
"""
psycho_simLR(T)

generate a psychometric curve based on simulated data returned from do_sim (or something similar) 
"""
psycho_simLR(T; β=0.3, x₀=NaN, as_choice=false) = begin
	lotteries, control, muscimol = T
	control = control[1,:,:] - control[2,:,:] # lottery - surebet
	muscimol = muscimol[1,:,:] - muscimol[2,:,:] # lottery - surebet
	
	
	if isnan(x₀)
		x₀ = quantile(control[:],0.5)
	end
	softmax(x) = 1 / (1 + exp(-β * (x - x₀)))
	if as_choice
		choose(x) = (rand() < 1 / (1 + exp(-β * (x - x₀))))
	else
		choose = softmax
	end
	X = reshape(repeat(lotteries, inner=size(control,1)), size(control,1),6)
	Y_control = choose.(control)
	Y_musc = choose.(muscimol)
	#lm = glm(@formula(chose_lottery ~ value), DataFrame(chose_lottery=Y_control[:], value=X[:]),Bernoulli)

	f1(p) = fitlogit(p, X[:] .+ 0.0, Y_control[:])
	o1 = optimize(f1, [51.0 ,0.1])
	p1 = o1.minimizer

	f2(p) = fitlogit(p, X[:] .+ 0.0, Y_musc[:])
	o2 = optimize(f2, [51.0 ,0.1])
	p2 = o2.minimizer
	
	
	return (X, Y_control, Y_musc, p1, p2)
	
	
	
end
	
		
		

# ╔═╡ 4a05ea33-5da0-4593-884d-51c26eda74f0
"""
do_sim2
"""
do_sim2(;boots=20, opts=OptionsLR(), W=[0.5 -1; -1 0.5], inh=([1],0.3), seed=1234) = begin
	Random.seed!(seed)
	values = [12*x for x in [0,1,2,4,8,16]]
	cntrl = fill(0.0, 2 , boots,  length(values))
	opto = similar(cntrl)
	for (vi,value) = enumerate(values)
			C = hcat(ThreadsX.map(x->simM2trace(([1],1.0); opts,W, value)[3][end,:], 1:boots)...)
			cntrl[:, :,vi] = C
			C = hcat(ThreadsX.map(x->simM2trace(inh; value, opts,W)[3][end,:],  1:boots)...)
			opto[:, :, vi] = C
	end
	return (values, cntrl, opto)
end

# ╔═╡ c3a82108-2219-4f1c-9884-117919a80af4
simLR_out = let
	# opto = ([1], 1.0)
	# T = simM2trace(opto; value=50)
	# p1 = plot(T[1], T[3], xlim=(-1,1.5))
	TT1 = do_sim2(W=WLR, opts=olr, inh=([1],0.2), boots=50)
	TT2 = do_sim2(W=WLR, opts=olr, inh=([1,2],0.2), boots=50)
	(TT1, TT2)
end; md"#### Run unilateral and bilateral inactivations"

# ╔═╡ 31ce1452-7669-44bd-a4eb-97d35f2421eb
begin
 	v,c,m = simLR_out[1]
	#[ci(c) for c in eachcol(c[1,:,:] - c[2,:,:])]
end

# ╔═╡ d5597ed9-8324-4ec2-994d-5cecc9412f60
quantile.(eachcol(c[1,:,:] - c[2,:,:]),[0.025 0.975])

# ╔═╡ b880b07b-bf5a-4cde-9afb-7e6f98d15932
ooo=let
	vals, cntrl, inh = simLR_out[1]
	q = (cntrl, inh)
	mean(q[2][2,:,:], dims=1), vals
	mean(cntrl[1,:,:] .- cntrl[2,:,:], dims=1)
end

# ╔═╡ bab5314c-0f06-4bba-88ee-8f320f44c6ec
let
	TT1, TT2 = simLR_out
	v,c,m = TT1
	mean(c[1,:,:], dims=1)
end

# ╔═╡ 52ac3898-6360-4d6c-b4ab-21095ce208b3
simM2trace(([1],1.0),value=20, opts=olr,W = WLR )[3][end,:]

# ╔═╡ 73052deb-5977-47eb-82aa-5c30afcd8715
foo(;a=10, b=10, c=10) = a + b + c

# ╔═╡ b900c1b9-b60e-423e-940e-434efc6009ea
let
	a = [[1 2 ],[3 4]]
	vcat(a...)
end

# ╔═╡ 7f17f7fc-5b27-4a10-9075-af9ad33633f8
cntrl_clr() = colorant"rgb(131,139,139)"

# ╔═╡ d132f0c0-fc0e-44ac-abe3-6ddfaf383bf1
fof_clr() = colorant"#a740ff"

# ╔═╡ cbed5b06-d59c-47a3-8f7f-340de9be56ee
"""
plot_sim(T)

T is the 2nd output of `do_sim`.

This function plots a summary of the inactivation results along with a power law fit to the relationship between U and V.


"""
plot_sim(T) = begin
	
	vals, cntrl, fof = T
	# Plot the U=V^\rho curves. Would be better to call the fit function locally or to pass it in instead of globals. :shrug:
	p1 = o1.minimizer
	p2 = o2.minimizer
	p = plot(x->x^p1[1] + p1[2], color=cntrl_clr(),label="",xlim=(0,200))
	plot!(x->x^p2[1] + p2[2], color=fof_clr(), label="")
	
	# Plots the error bars for cntrl
	scatter!(vals,mean(eachrow(cntrl)), yerr=ci(cntrl), lc=cntrl_clr(),
		color=cntrl_clr(), ms=0,label="", legend=:topleft)
	
	# Plots the markers for cntrl
	scatter!(vals,mean(eachrow(cntrl)), 
		color=cntrl_clr(), label="Control")
	
	# Plot the error bars for muscimol
	jitter = 3
	scatter!(vals.+jitter,mean(eachrow(fof)), yerr=ci(fof), lc=fof_clr(),msc=:auto,
		color=fof_clr(),ms=0,label="", legend=:topleft,
		xlab="Lottery Magnitude\n(μL of water)",
		ylab="Firing Rate (Hz)",
		leftmargin=40px,
		bottommargin=40px)
	
	scatter!(vals.+jitter,mean(eachrow(fof)), 
		color=fof_clr(), label="Muscimol", shape=:diamond)

	hline!([sb_val],ls=:dot,color=:blue, label="")
	
	p
end

# ╔═╡ eea4905d-9fa6-46ed-8b13-dfa426cc6106
"""
Plot psychometrics based on simulation
"""
plot_psycho(T) = begin
    v,c,m,pc,pm=psycho_sim(T)
	p=scatter(v[1,:], mean(eachrow(c)), yerr=ci(c), color=cntrl_clr(), ms=0)
	scatter!(v[1,:], mean(eachrow(c)), color=cntrl_clr())
	scatter!(v[1,:].+3, mean(eachrow(m)), yerr=ci(m), color=fof_clr(), ms=0)
	scatter!(v[1,:].+3, mean(eachrow(m)), color=fof_clr(), 
		shape=:diamond,
		legend=false)
	fitx = 0:1:200
	plot!(fitx, logit(pc, fitx), color=cntrl_clr())
	plot!(fitx, logit(pm, fitx), color=fof_clr(),
		xlabel="Lottery Magnitude",
		ylabel="P(Choose Lottery)",
		leftmargin=40px,
		bottommargin=40px,
		yticks=[0:0.5:1;], ylim=(-0.1,1.1)
		
	)
	
	p
end

# ╔═╡ 9252406d-b27d-420c-8029-a68031078910
"""
plot_leaky(TL)

plot the results of the 
"""
plot_leaky(TL) = begin
	t = TL[1]
	D = map(x->mean(eachrow(x)), TL[2:end])
	ls = :dash
	# cntrl_low, cntrl_high, opto_low, opto_high)
	i=1
	
	plot(t, D[i]; label="C ↑" ,color=cntrl_clr())
	scatter!(t, D[i]; yerr=ci(TL[i+1]), lc=cntrl_clr(),
		color=cntrl_clr(), ms=0,label="")

	i+=1
	plot!(t.+i, D[i];ls, label="C ↓" ,color=cntrl_clr())
	scatter!(t.+i, D[i]; yerr=ci(TL[i+1]), lc=cntrl_clr(),
		color=cntrl_clr(), ms=0,label="")
	
	i+=1
	color=fof_clr()
	plot!(t.+i, D[i]; label="O ↑" ,color)
	scatter!(t.+i, D[i]; yerr=ci(TL[i+1]), lc=color,
		color, ms=0,label="")

	i+=1
	plot!(t.+i, D[i];ls, label="O ↓" ,color)
	scatter!(t.+i, D[i]; yerr=ci(TL[i+1]), lc=color,
		color, ms=0,label="")
	

end

# ╔═╡ e1b770ec-f40d-485c-bafa-b09e8ccad7c6
trace_plot(;value=12*8, seed=1234) = begin
	Random.seed!(seed)
	t,μ,U,W = simMint(0;value)
	α = fill(0.1, (1,size(U,2)))
	α[1] = .5
	p = plot(t,U, label="", legend=:topleft, alpha=α,
	xlab="Time from Sound onset (s)", color=cntrl_clr(),
	ylab="Firing Rate (Hz)", size=(500,500), xlim=(-0.25,1.0),
)
	tf,μf,Uf = simMint(0.2; value,W)
	
	
	plot!(p,tf, Uf, label="", color=fof_clr(), alpha=0.1)
	
	plot!(t, mean(eachcol(U)), color=cntrl_clr(), linewidth=3, label="Control")
	plot!(t, mean(eachcol(Uf)), color=fof_clr(), linewidth=3, label="Muscimol")
	return p
end; md"#### trace_plot"

# ╔═╡ bda5109e-7f98-45ec-a647-f6a248c376ea
begin
p2 =plot(
		plot(trace_plot(seed=my_seed),xticks=[0,0.5,1]), 
		plot_sim(T), 
		plot_psycho(T),
		foreground_color_legend = nothing, 
		layout=grid(1,3), size=(600,220),
		guidefontsize=10, leftmargin=40px,
		rightmargin=10px
)
md"""
$(p2)

**Left: Single Trial Simulation.** The thin lines show the activity on a single trial of the 6 nodes in the network. The thick lines show the average activity of the network in control (gray) and FOF silencing (purple) conditions. The slightly darker gray line is the FOF activity during the control condition.

**Middle: Effect of muscimol.** The dots with error bars represent the network responses across 20 trials (per lottery). In gray are the reponses on control sessions, the purple points are data from when the FOF node is silenced. The lines represent power law fits of activity to lottery value, _Hz = uLᵖ + κ_. The dotted blue line represents the remembered value of the surebet. When activity is above that line, the network is more likely to choose the lottery, when activity is below that line, the network is more likely to choose the sure bet. 

**Right: Simulated psychometric curves.** We applied a logistic choice function to the activity at the time of the go cue (the same activity plotted in the middle panel). The network activity in the middle panel results in psychometric curves that qualitatively match the data from muscimol and optogenetic silencing.

"""
end

# ╔═╡ 8c48ffe4-f8c2-4338-8dc3-0d4349e4469c
begin
	print_it
	
	savefig(plot(p2, size=(750,250)),ENV["HOME"]*"/papers/risk-musc/overleaf-risk-musc/figs/bio_6_$(my_seed).pdf")
end

# ╔═╡ 1d26b5dc-e79f-4289-8a67-6166988c8def
# A way to keep track of some reasonable seeds
manyp = plot(map(x->trace_plot(seed=x),[709,588,105,56,125,131,205,42,637])..., legend=:topleft, size=(800,500)); md"""
$(manyp)
Some other simluations of the same network (different random weight matrix.)

"""

# ╔═╡ a09c0ee7-e360-4ff5-859f-c5987475e17e
"""
`plot_simLR(TT)`

Where TT is a 3-tule with (vals, cntrl, inhibition) 
"""
plot_simLR(TT,str) = begin
	μ(x) = vec(mean(x, dims=1))
	c(x) = [0.6 .* x, x]
	m() = [:square :circle]
	# @show TT[1] d(TT[2])
	vals, cntrl, inh = TT

	p = plot(;
		xlabel="Lottery Magnitude",
		ylabel="Hz",
		foreground_color_legend = nothing,
		legend=:topleft
	)
	q = (cntrl, inh)
	clrs = (cntrl_clr(), fof_clr())
	jitter = (0,3)
	labels = (["cntl-L" "cntl-SB"], str=="lottery" ? ["inh-L" "cntl-SB"] : ["inh-L" "inh-SB"])
	for i ∈ [1,2] # i indexes cntrl, inhibition
		for j ∈ [1,2] # j indexes lottery vs. surebet node
			y = q[i][j,:,:]
			total_jitter = jitter[i] .+ j*4 
			clr = clrs[i]
			#j == 2 && (clr *= 1.3)
			scatter!(vals .+ total_jitter, μ(y), yerr=ci(y), lc=clr,
		 		ms=0,label="")
			scatter!(vals .+ total_jitter, μ(y), 
				color=clr,
				marker=m()[j], 
				label=labels[i][j])
			plot!(vals .+ total_jitter, μ(y), 
				color=clr,
				label="")
			
			
		end
	end
	p
end

# ╔═╡ c5c98f83-71de-4eb5-8edd-f867b73901d5
"""
plot_psychoLR(T)

Plot psychometrics based on simulation for 2-node LR model
"""
plot_psychoLR(T) = begin
    v,c,m,pc,pm=psycho_simLR(T)
	p=scatter(v[1,:], mean(eachrow(c)), yerr=ci(c), color=cntrl_clr(), ms=0)
	scatter!(v[1,:], mean(eachrow(c)), color=cntrl_clr())
	scatter!(v[1,:].+3, mean(eachrow(m)), yerr=ci(m), color=fof_clr(), ms=0)
	scatter!(v[1,:].+3, mean(eachrow(m)), color=fof_clr(), 
		shape=:diamond,
		legend=false)
	fitx = 0:1:200
	plot!(fitx, logit(pc, fitx), color=cntrl_clr())
	plot!(fitx, logit(pm, fitx), color=fof_clr(),
		xlabel="Lottery Magnitude",
		ylabel="P(Choose Lottery)",
		leftmargin=40px,
		bottommargin=40px,
		yticks=[0:0.5:1;], ylims=(-0.1,1.1)
		
	)
	
	p
end

# ╔═╡ a7cc2e68-d2be-4623-bf82-ff0ce94a5e64
pLR = let
	TT1, TT2 = simLR_out
	TT = [x for x in TT2]
	TT[2] = TT1[2]
	p1=plot_simLR(TT1, "lottery")
	p2=plot_simLR(TT, "bilateral")
	
	pp1=plot_psychoLR(TT1)
	pp2=plot_psychoLR(TT)
	
	plot(p1,p2,pp1,pp2, layout=grid(2,2), size=(600,400); bottom_margin=10px)
end; md"""
$(pLR)
**Simulation of a two-node network, where the nodes represent the left and right FOF**. 

**Left Column: Unilateral Inhibition of the FOF contralateral to the lottery choice** Top: Control simulation (without inhibition) is plotted in gray. Simulation with inhibition of the lottery node in purple. The square points represent the mean (± 95% CI of the mean) response of the lottery node over $(size(simLR_out[1][2],2)) trials. The circles represent the responses of the surebet nodes. Bottom: Behavioral performance of the networks simulated in the top row. Grey circles represent a control simulation (with no perturbations) and Purple diamonds represent the choices for a unilateral inactivation of the side of the FOF contralateral to the lottery choice. This perturbation results in a vertical scaling of the psychometric curve as in Erlich 2011.

**Right Column: Bilateral Inhibition of the FOF** Similar to left column, but for bilateral inactivations. Not that for very bad lotteries, bilateral silencing _increases_ the number of lottery choices. Overall, bilateral silencing of this network decreases the slope of the psychometric curve. 
"""

# ╔═╡ 2ecccee9-08bc-439d-bb92-52229ac0349f
@bind fclr PlutoUI.ColorStringPicker()

# ╔═╡ d659c4d1-db14-4de4-b26a-b47bcedadc27
PlutoUI.TableOfContents(depth=6, include_definitions=true)

# ╔═╡ dc456191-8e25-4977-82fe-df270846926b
cntrl_clr(), fof_clr()

# ╔═╡ fe4b331c-a91a-441e-bce3-5a229d3597ae
let
	a = [1 2; 3  4]
	a[2,2] *= 2
	a
end

# ╔═╡ 0feed914-ca8d-43e5-9822-f705134e19f4
foo(x::Integer) = foo(x * 1.0)

# ╔═╡ 6b13fbdb-935e-4901-912d-4fa572b35ec1
foo(x::Number) = x+2

# ╔═╡ 8c4c5e87-d3c1-46ad-a540-4563e756322b
foo(x) = x - 1

# ╔═╡ 402fda2f-fe1e-4e7e-a0cd-e3234dc56fc2
foo(x::String) = "your number is $(x)"

# ╔═╡ 29bc97eb-2c91-4ad4-83da-3ec357e1d32a
foo(1)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
GLM = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
ThreadsX = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
UnPack = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"

[compat]
DataFrames = "~1.4.4"
GLM = "~1.8.1"
Optim = "~1.7.4"
Plots = "~1.38.2"
PlutoUI = "~0.7.49"
StatsBase = "~0.33.21"
ThreadsX = "~0.1.11"
UnPack = "~1.0.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.1"
manifest_format = "2.0"
project_hash = "1b114c63e01b7ac2635784d85fb86bfabb6d4df7"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "14c3f84a763848906ac681f94cf469a851601d92"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.28"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "7fe6d92c4f281cf4ca6f2fba0ce7b299742da7ca"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.37"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "00a2cccc7f098ff3b66806862d275ca3db9e6e5a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.5.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d4f69885afa5e6149d0cab3818491565cf41446d"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.4.4"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "c5b6685d53f933c11404a3ae9822afe30d522494"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.12.2"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "74911ad88921455c6afcad1eefa12bd7b1724631"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.80"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "9a0472ec2f5409db243160a8b030f94c380167a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.6"

[[deps.FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "04ed1f0029b6b3af88343e439b995141cb0d0b8d"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.17.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "a69dd6db8a809f78846ff259298678f0d6212180"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.34"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GLM]]
deps = ["Distributions", "LinearAlgebra", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "StatsModels"]
git-tree-sha1 = "884477b9886a52a84378275737e2823a5c98e349"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.8.1"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "387d2b8b3ca57b791633f0993b31d8cb43ea3292"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.71.3"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "5982b5e20f97bff955e9a2343a14da96a746cd8c"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.71.3+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "eb5aa5e3b500e191763d35198f859e4b40fff4a6"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.7.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.InvertedIndices]]
git-tree-sha1 = "82aec7a3dd64f4d9584659dc0b62ef7db2ef3e19"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.2.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "946607f84feb96220f480e0422d3484c49c00239"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.19"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "4d5917a26ca33c66c8e5ca3247bd163624d35493"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "6503b77492fd7fcb9379bf73cd31035670e3c509"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "1903afc76b7d01719d9c30d3c7d501b61db96721"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.4"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "8175fc2b118a3755113c8e68084dc1a9e63c61ee"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.3"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "5b7690dd212e026bbab1860016a6601cb077ab66"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.2"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "a99bbd3664bb12a775cda2eba7f3b2facf3dad94"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eadad7b14cf046de6eb41f13c9275e5aa2711ab6"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.49"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "96f6db03ab535bdb901300f88335257b0018689d"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "97aa253e65b784fd13e83774cadc95b38011d734"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.6.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase", "SnoopPrecompile"]
git-tree-sha1 = "e974477be88cb5e3040009f3767611bc6357846f"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.11"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "e681d3bfa49cd46c3c161505caddf20f0e62aaa9"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.ShiftedArrays]]
git-tree-sha1 = "503688b59397b3307443af35cd953a13e8005c16"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "2.0.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "6954a456979f23d05085727adb17c4551c19ecd1"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.12"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "ab6083f09b3e617e34a956b43e9d51b824206932"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.1.1"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "a5e15f27abd2692ccb61a99e0854dfb7d48017db"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.6.33"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadsX]]
deps = ["ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "Setfield", "SplittablesBase", "Transducers"]
git-tree-sha1 = "34e6bcf36b9ed5d56489600cf9f3c16843fa2aa2"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.11"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "c42fa452a60f022e9e087823b47e5a5f8adc53d5"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.75"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─74fb6a63-48a4-4f52-bf49-dfcc70785f66
# ╟─0a41f9e0-36d8-406e-8fca-f6f277b03d1c
# ╠═0d028ece-4112-4518-b320-f896704c97ab
# ╠═bda5109e-7f98-45ec-a647-f6a248c376ea
# ╠═1d26b5dc-e79f-4289-8a67-6166988c8def
# ╠═2263c34f-1481-40f0-9893-32334ce65009
# ╟─a4a67fe2-3dcd-4d0f-a514-d15ab736e01e
# ╟─8c48ffe4-f8c2-4338-8dc3-0d4349e4469c
# ╟─02ac31e6-e9d1-4946-bbd1-9ed802310d63
# ╠═a7cc2e68-d2be-4623-bf82-ff0ce94a5e64
# ╠═eec8c196-9808-44db-a65b-dc3c52b8c04e
# ╠═31ce1452-7669-44bd-a4eb-97d35f2421eb
# ╠═d5597ed9-8324-4ec2-994d-5cecc9412f60
# ╟─d38ef57a-2519-4f95-b5a2-8d890b29915e
# ╠═f59fce3a-3273-11ec-3fe8-9fed7d8c159e
# ╠═c28c212e-3f37-4ad5-ae4b-9377784a22f0
# ╠═40f6aa19-583b-4e88-93e7-fabf16fd082f
# ╟─f4ba91f5-aba0-46a3-a8b6-b2da0fb4bbea
# ╠═cbed5b06-d59c-47a3-8f7f-340de9be56ee
# ╠═eea4905d-9fa6-46ed-8b13-dfa426cc6106
# ╟─917f1cd3-6545-42de-a2ea-e48b9d6f9665
# ╠═3610a639-8a8c-4af6-b161-0c5badfc3e8e
# ╟─9252406d-b27d-420c-8029-a68031078910
# ╟─903f6ce1-66e8-470b-815f-2a6a69581d9a
# ╠═8d03feb6-8d94-4765-914a-266b7b04034e
# ╟─f0e68369-c342-43e1-80d3-71ef27eede36
# ╟─10306590-a806-4708-a823-07d11d5c1c2d
# ╠═3aee647e-6150-448f-af48-6ddcadd4f688
# ╟─5789c252-9ffa-479c-8bd1-bad310ee9005
# ╠═752ea43e-5f42-410c-acfe-2dcc795b4c89
# ╠═e1b770ec-f40d-485c-bafa-b09e8ccad7c6
# ╠═030cea43-fecf-417a-9e1d-7a6970569234
# ╠═923cacc6-f84d-44c3-829b-9b82d121a80c
# ╟─aaebd031-a8de-40d8-a2da-8c9408d89e0d
# ╠═65364d78-2b7e-4fc6-80e8-d4a1af740676
# ╠═8eb705a5-f23e-4c60-bfff-5fd2752b363e
# ╠═1d488acb-f9d3-4de8-b193-2016159e8ea2
# ╠═c3a82108-2219-4f1c-9884-117919a80af4
# ╠═c28e9c6b-5ffa-40f5-be7b-642bca8d6c65
# ╠═b880b07b-bf5a-4cde-9afb-7e6f98d15932
# ╠═90b5837b-c763-47ac-91a8-1354c3551e2d
# ╟─a09c0ee7-e360-4ff5-859f-c5987475e17e
# ╠═c5c98f83-71de-4eb5-8edd-f867b73901d5
# ╠═05e638d3-6b71-44f4-a87d-d920bee04b1b
# ╠═bab5314c-0f06-4bba-88ee-8f320f44c6ec
# ╠═4a05ea33-5da0-4593-884d-51c26eda74f0
# ╠═52ac3898-6360-4d6c-b4ab-21095ce208b3
# ╠═73052deb-5977-47eb-82aa-5c30afcd8715
# ╠═b900c1b9-b60e-423e-940e-434efc6009ea
# ╠═7f17f7fc-5b27-4a10-9075-af9ad33633f8
# ╠═d132f0c0-fc0e-44ac-abe3-6ddfaf383bf1
# ╠═2ecccee9-08bc-439d-bb92-52229ac0349f
# ╠═d659c4d1-db14-4de4-b26a-b47bcedadc27
# ╠═dc456191-8e25-4977-82fe-df270846926b
# ╠═fe4b331c-a91a-441e-bce3-5a229d3597ae
# ╠═0feed914-ca8d-43e5-9822-f705134e19f4
# ╠═6b13fbdb-935e-4901-912d-4fa572b35ec1
# ╠═8c4c5e87-d3c1-46ad-a540-4563e756322b
# ╠═402fda2f-fe1e-4e7e-a0cd-e3234dc56fc2
# ╠═29bc97eb-2c91-4ad4-83da-3ec357e1d32a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
