using DelimitedFiles
using ShiftedArrays
using Distributed

# to use in parallel settings
# addprocs(1)

@everywhere begin
    Base.include(Main, "multipop_MANTIS.jl")
    using ShiftedArrays
    function count_unique_local_minima(vector, digits=nothing)
        # -1 = extinct, 0 = stable, 1 = unconverged, >1 = chaos/cycles
        vector = isnothing(digits) ? vector : round.(vector, digits=digits)
        # if the last element is 0, the strain is extinct, regardless of other dynamics
        if vector[end] == 0 return(-1) end
        # pick out the local minima (use isless bc handles missing values introduced by lead/lag
        # silently: missing > all numbers)
        minima = vector[isless.(vector, ShiftedArrays.lag(vector)) .& isless.(vector, ShiftedArrays.lead(vector))]
        number_unique_minima = length(unique(minima))
        # if there are no minima, the strain is constant or has a period proportional to the stepsize (unlikely)
        # if only one minima, the vector is monotonic (implies failure to converge)
        if number_unique_minima <= 1 return(number_unique_minima) end
        # now check for monotonicity (read: failure to converge) in cyclical dynamics (i.e. length(minima) > 1)
        if all(((minima - ShiftedArrays.lag(minima)) .<= 0)[2:end]) |
             all(((minima - ShiftedArrays.lag(minima)) .>= 0)[2:end])
            return(1)
        end
        # if not caught by any of the earlier conditions, then the dynamics are chaotic/cyclical
        return(number_unique_minima)
    end
    # for creating state-space figures
    function clean_sse_timeseries_output(tmp, params, nstrains)
        ntimepoints = size(tmp, 1)
        tmp = DataFrames.stack(tmp, 1:size(tmp, 2) - 1, size(tmp, 2))
        tmp.variable = vcat([[x for y in 1:ntimepoints] for x in 1:nstrains]...)
        # note we round here to improve correct classification of local extrema
        tmp = DataFrames.by(tmp, [:variable], :value => x -> count_unique_local_minima(x, 8))
        tmp[!, :gamma] .= params[1]
        tmp[!, :sigma] .= params[2]
        tmp[!, :R0] .= params[3]
        tmp[!, :mu] .= params[4]
        return(tmp)
    end
    function clean_vmr_timeseries_output(tmp, params, npopulations, nstrains)
        ntimepoints = size(tmp, 1)
        tmp = DataFrames.stack(tmp, 1:size(tmp, 2) - 1, size(tmp, 2))
        tmp.variable = vcat([[x for y in 1:2*ntimepoints] for x in 1:nstrains]...)
        tmp.population = vcat([vcat(fill(1, ntimepoints), fill(2, ntimepoints)) for x in 1:nstrains]...)
        # note we round here to improve correct classification of local extrema
        tmp = DataFrames.by(tmp, [:variable, :population], :value => x -> count_unique_local_minima(x, 8))
        tmp.chi = fill(params[1], nstrains * npopulations)
        tmp.gamma = fill(params[2], nstrains * npopulations)
        tmp.sigma = fill(params[3], nstrains * npopulations)
        tmp.R0 = vcat(fill(params[4:5], nstrains)...)
        return(tmp)
    end
end

function state_space_exploration(;strainstruct = [2 2],
                                 gamma         = [0.44 0.55 0.66 0.77 0.88]
                                 sigma         = [1 2 4 8 16 32 64]
                                 rnaught       = range(1, length=40, stop=5)
                                 mu            = range(0.01, length=40, stop=0.25)
                                 maxtime       = 10000.0)
    params = [[γ, σ, R0, μ] for γ in gamma for σ in sigma for R0 in rnaught for μ in mu]
    nstrains = trunc(Int, prod(strainstruct))
    results = pmap(x -> clean_eps1p_timeseries_output(
                            runMANTIS(strainstructure=strainstruct,
                                      tmax=maxtime,
                                      tstep=trunc(Int, 0.8*maxtime):maxtime,
                                      beta=x[3] * x[2], gamma=x[1], sigma=x[2],
                                      chi=0.0, mu=x[4])["timeseries"][:,vcat((nstrains + 1):(2 * nstrains), end)], x, nstrains),
                   params)
    return(convert(Matrix, vcat(results...)))
end

# provide parameter combinations as a vector of tuples: (γ, σ, R01, R02)
function variable_movement_rate(;param_combinations = [(0.55, 32, 2, 5),
                                                       (0.55, 32, 5, 2),
                                                       (0.66,  8, 2, 5),
                                                       (0.66,  8, 5, 2),
                                                       (0.77,  4, 3, 5),
                                                       (0.77,  4, 5, 3)],
                                strainstruct        = [2 2],
                                chi                 = 10.0 .^ -(1:0.5:7),
                                maxtime             = 10000.0,
                                mu                  = 0.15)

    params = [[χ param_comb...] for χ in chi for param_comb in param_combinations]
    nstrains = trunc(Int, prod(strainstruct))
    npopulations = 2
    results = pmap(x -> clean_vmr_timeseries_output(
                            runMANTIS(strainstructure=strainstruct,
                                      tmax=maxtime,
                                      tstep=trunc(Int, 0.8*maxtime):maxtime,
                                      beta=x[4:5] * (x[3] .+ mu), gamma=x[2], sigma=x[3],
                                      chi=[-x[1] x[1]; 0 0], mu=mu .- [x[1]; 0])["timeseries"][:,vcat((npopulations * nstrains + 1):(2 * npopulations * nstrains), end)], x, npopulations, nstrains),
                   params)
    return(convert(Matrix, vcat(results...)))
end
