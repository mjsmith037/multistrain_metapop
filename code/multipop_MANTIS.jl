using DataFrames
using Random
using StaticArrays
using DifferentialEquations

# legacy seasonality implementation from MANTIS framework
seasonality(t, β, ε) = (1 + ε * sin(π * t)^6) .* β

# governing equations
dydt(y, z, w, λ, γ, σ, μ, χ) = λ .* ((1 .- w) + (1 .- γ) .* (w - z)) - σ .* y - μ .* y + χ' * y
dzdt(z, λ, μ, χ) = λ .* (1 .- z) - μ .* z + χ' * z
dwdt(w, ssλ, μ, χ) = ssλ .* (1 .- w) - μ .* w + χ' * w

function mantis!(deriv, state, params, t)
     y = @view state[:,:,1];  z = @view state[:,:,2];  w = @view state[:,:,3]
    dy = @view deriv[:,:,1]; dz = @view deriv[:,:,2]; dw = @view deriv[:,:,3]
    # adjust β to address seasonality
    βi = seasonality(t, params.β, params.ε)
    # force of infection in y,z equations
    λ = βi .* y
    # force of infection for strains that share epitopes with each strain
    # [i.e. force of infection for each virus in w equation]
    ssλ = hcat([sum(βi[:,x] .* view(y, :, params.sharedalleles[:,x]), dims=2) for x in 1:size(y, 2)]...)
    # calculate the derivatives
    dy .= dydt(y, z, w, λ, params.γ, params.σ, params.μ, params.χ)
    dz .= dzdt(z, λ, params.μ, params.χ)
    dw .= dwdt(w, ssλ, params.μ, params.χ)
end

# small macro to get a variable name as a string
macro Name(arg)
   string(arg)
end

# resize parameters by replicating across populations, strains, or both
function parametersizecheck(parameter, parametername::String, npops::Int64, nstrains::Int64)
    if length(parameter) == 1
        resized_parameter = SMatrix{npops, nstrains}(fill(parameter, npops, nstrains))
    # note, if npops == nstrains, this formulation assumes differences
    # are between populations rather than between strains
    elseif length(parameter) == npops
        resized_parameter = SMatrix{npops, nstrains}(hcat(fill(reshape(parameter, npops, 1), nstrains)...))
    elseif length(parameter) == nstrains
        resized_parameter = SMatrix{npops, nstrains}(vcat(fill(reshape(parameter, 1, nstrains), npops)...))
    elseif length(parameter) == npops * nstrains
        resized_parameter = SMatrix{npops, nstrains}(reshape(parameter, npops, nstrains))
    else error("Invalid $parametername value") end
    return(resized_parameter)
end

function runMANTIS(;strainstructure=[2 2], tmax=1000, tstep=1:1000,
                   beta=16, gamma=0.66, epsilon=0, sigma=8, mu=0.1,
                   chi=[-0.05 0.0; 0.0 -0.05], rseed=1, initvals::Any="random")

    # if movement is provided as a single value (i.e. movement
    # of a single-population system) make it a matrix
    chi = isempty(size(chi)) ? SMatrix{1,1}(chi) : SMatrix{size(chi)...}(chi)

    # enumerate the possible strains based on the number of loci, alleles
    strains = reshape(collect(Iterators.product([1:x for x in strainstructure]...)), :, 1)

    npops = size(chi, 1)
    nstrains = length(strains)

    # generate initial conditions if necessary (can also be provided as
    # Array{npops, nstrains} corresponding to initial y values)
    if typeof(initvals) == String && initvals == "random"
        # Random.seed!(rseed)
        initvals = rand(npops, nstrains) |> x->x./10 #sum(x, dims=2)
    end

    # which strains share alleles?
    sharedalleles = hcat([[any(x .== y) for x in strains] for y in strains]...)

    # reshape the initial conditions and adjust w initial conditions based on shared alleles
    initcond = cat(initvals, initvals,
                   min.(hcat([sum(initvals[:,sharedalleles[:,x]], dims=2) for x in 1:nstrains]...), 1),
                   dims=3)

    # reshape all parameters prior to analysis and collect into a named tuple
    params = (β = parametersizecheck(beta, @Name(beta), npops, nstrains),
              σ = parametersizecheck(sigma, @Name(sigma), npops, nstrains),
              μ = parametersizecheck(mu, @Name(mu), npops, nstrains),
              γ = parametersizecheck(gamma, @Name(gamma), npops, nstrains),
              ε = epsilon, χ = chi, sharedalleles = sharedalleles)

    # run the simulation and return a timeseries as well as the parametrization
    prob = ODEProblem(mantis!, initcond, (0.0, tmax), params)
    sol = solve(prob, saveat=tstep, abstol=1e-12, reltol=1e-12,
                isoutofdomain=(u,p,t) -> any(x -> (x < 0) | (x > 1), u))
    return(Dict("timeseries" => filter(row -> row[end] ∈ tstep,
                                       DataFrame(hcat([[reshape(u, (1,prod(size(u))))..., t]
                                                           for (u,t) in tuples(sol)]...)')),
                "parameters" => params))
end
