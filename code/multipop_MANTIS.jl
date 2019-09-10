import DataFrames
import Random
import StaticArrays
using DifferentialEquations

seasonality(t, β, ε) = (1 + ε * sin(π * t)^6) .* β

dydt(y, z, w, λ, γ, σ, μ, χ) = λ .* ((1 .- w) .+ (1 .- γ) .* (w .- z)) - σ .* y .- μ .* y .+ χ' * y
dzdt(z, λ, μ, χ) = λ .* (1 .- z) .- μ .* z .+ χ' * z
dwdt(w, ssλ, μ, χ) = ssλ .* (1 .- w) .- μ .* w .+ χ' * w

function mantis!(deriv, state, params, t)
     y = @view state[:,:,1];  z = @view state[:,:,2];  w = @view state[:,:,3]
    dy = @view deriv[:,:,1]; dz = @view deriv[:,:,2]; dw = @view deriv[:,:,3]
    # adjust β to address seasonality
    βi = seasonality(t, params.β, params.ε)
    # force of infection in y,z equations
    # λ = βi .* y
    λ = βi .* y
    # force of infection for strains that share epitopes with each strain
    # [i.e. force of infection for each virus in w equation]
    ssλ = hcat([sum(βi[x] .* view(state, :, params.sharedalleles[x,:], 1), dims=2) for x in 1:size(state, 2)]...)
    # calculate the derivatives
    dy .= dydt(y, z, w, λ, params.γ, params.σ, params.μ, params.χ)
    dz .= dzdt(z, λ, params.μ, params.χ)
    dw .= dwdt(w, ssλ, params.μ, params.χ)
end

function runMANTIS(;strainstructure=[2 2], tmax=1000.0, tstep=missing,
                   beta=4.0, gamma=0.72, epsilon=0.0, sigma=1.0, mu=0.0,
                   chi=[-0.02 0.0; 0.0 -0.02], rseed=1, initvals::Any="random")
    # strainstructure=[2 2];tmax=100.0;tstep=missing;beta=40.0;gamma=0.72;
    # epsilon=0.0;sigma=10.0;mu=0.0;chi=[-0.02 0.0; 0.0 -0.02];rseed=1;initvals="random"

    chi = isempty(size(chi)) ? StaticArrays.SMatrix{1,1}(chi) : StaticArrays.SMatrix{size(chi)...}(chi)

    strains = reshape(collect(Iterators.product([1:x for x in strainstructure]...)), :, 1)

    npops = size(chi, 1)
    nstrains = length(strains)

    if typeof(initvals) == String && initvals == "random"
        # Random.seed!(rseed)
        Random.seed!(trunc(Int, rand() * 1000))
        initvals = rand(npops, nstrains) |> x->x./sum(x, dims=2)
    end

    # initcond = StaticArrays.Size((size(initvals)..., 3))([initvals..., initvals..., min.(3*initvals, 1)...])
    initcond = cat(initvals, initvals, min.(3*initvals, 1), dims=3)

    if length(beta) == 1
        beta_resized = fill(beta, npops, nstrains)
    elseif length(beta) == npops
        beta_resized = hcat(fill(reshape(beta, npops, 1), nstrains)...)
    elseif length(beta) == nstrains
        beta_resized = vcat(fill(reshape(beta, 1, nstrains), npops)...)
    elseif length(beta) == npops * nstrains
        beta_resized = reshape(beta, npops, nstrains)
    else error("Invalid beta value") end

    if length(sigma) == 1
        sigma_resized = fill(sigma, npops, nstrains)
    elseif length(sigma) == npops
        sigma_resized = hcat(fill(reshape(sigma, npops, 1), nstrains)...)
    elseif length(sigma) == nstrains
        sigma_resized = vcat(fill(reshape(sigma, 1, nstrains), npops)...)
    elseif length(sigma) == npops * nstrains
        sigma_resized = reshape(sigma, npops, nstrains)
    else error("Invalid sigma value") end

    if length(mu) == 1
        mu_resized = StaticArrays.@SMatrix fill(mu, npops, nstrains)
    elseif length(mu) == npops
        mu_resized = hcat(fill(reshape(mu, npops, 1), nstrains)...)
    elseif length(mu) == nstrains
        mu_resized = vcat(fill(reshape(mu, 1, nstrains), npops)...)
    elseif length(mu) == npops * nstrains
        mu_resized = reshape(mu, npops, nstrains)
    else error("Invalid mu value") end

    if length(gamma) == 1
        gamma_resized = StaticArrays.@SMatrix fill(gamma, npops, nstrains)
    elseif length(gamma) == npops
        gamma_resized = hcat(fill(reshape(gamma, npops, 1), nstrains)...)
    elseif length(gamma) == nstrains
        gamma_resized = vcat(fill(reshape(gamma, 1, nstrains), npops)...)
    elseif length(gamma) == npops * nstrains
        gamma_resized = reshape(gamma, npops, nstrains)
    else error("Invalid gamma value") end

    params = (β = beta_resized,
              γ = gamma_resized,
              ε = epsilon,
              σ = sigma_resized,
              μ = mu_resized,
              χ = chi,
              sharedalleles = hcat([[any(x .== y) for x in strains] for y in strains]...))

    prob = ODEProblem(mantis!, initcond, (0.0, tmax), params)
    if ismissing(tstep)
        sol = solve(prob, save_everystep=false)
        return(sol)
    else
        sol = solve(prob, AutoTsit5(Rosenbrock23()), reltol=1e-12, abstol=1e-12, saveat=tstep)
        tuplestosave = tuples(sol)#[trunc(Int, tmax / tstep * 0.9):end]
        return(Dict("timeseries" => DataFrames.DataFrame(hcat([[reshape(u, (1,prod(size(u))))..., t] for (u,t) in tuplestosave]...)'),
                    "parameters" => params))
    end
end
