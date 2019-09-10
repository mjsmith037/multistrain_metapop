using Distributed

while nprocs() < length(Sys.cpu_info())
    addprocs(1)
end

@everywhere begin
    Base.include(Main, "multipop_MANTIS.jl")
    import ShiftedArrays
    function is_local_minima(vector)
        vector = unique(round.(vector, digits=8))
        isless.(vector, ShiftedArrays.lag(vector)) .& isless.(vector, ShiftedArrays.lead(vector))
    end
    function clean_eps_timeseries_output!(tmp, params, npopulations, nstrains)
        ntimepoints = size(tmp, 1)
        tmp = DataFrames.stack(tmp, 1:size(tmp, 2) - 1, size(tmp, 2))
        tmp.variable = vcat([[x for y in 1:2*ntimepoints] for x in 1:nstrains]...)
        tmp.population = vcat([vcat([1 for y in 1:ntimepoints], [2 for y in 1:ntimepoints]) for x in 1:nstrains]...)
        tmp = DataFrames.by(tmp, [:variable, :population], :value => x -> sum(unique(round.(x, digits=8)) .> 0))
        unique!(tmp)
        tmp = DataFrames.by(tmp, [:variable, :population], :value_function => x -> x[is_local_minima(x)])
        tmp.gamma = vcat([vcat(params[1:2]...) for x in 1:nstrains]...)
        tmp.sigma = vcat([vcat(params[3:4]...) for x in 1:nstrains]...)
        tmp.R0 = vcat([params[5] for x in 1:nstrains*npopulations]...)
        return(tmp)
    end
    function clean_eps1p_timeseries_output!(tmp, params, nstrains)
        ntimepoints = size(tmp, 1)
        tmp = DataFrames.stack(tmp, 1:size(tmp, 2) - 1, size(tmp, 2))
        tmp.variable = vcat([[x for y in 1:ntimepoints] for x in 1:nstrains]...)
        tmp = DataFrames.by(tmp, [:variable], :value => x -> sum(unique(round.(x, digits=8)) .> 0))
        unique!(tmp)
        tmp = DataFrames.by(tmp, [:variable], :value_function => x -> x[is_local_minima(x)])
        tmp[!, :gamma] .= params[1]
        tmp[!, :sigma] .= params[2]
        tmp[!, :R0] .= params[3]
        tmp[!, :mu] .= params[4]
        return(tmp)
    end
    function clean_vmr_timeseries_output!(tmp, params, npopulations, nstrains)
        ntimepoints = size(tmp, 1)
        tmp = DataFrames.stack(tmp, 1:size(tmp, 2) - 1, size(tmp, 2))
        tmp.variable = vcat([[x for y in 1:2*ntimepoints] for x in 1:nstrains]...)
        tmp.population = vcat([vcat([1 for y in 1:ntimepoints], [2 for y in 1:ntimepoints]) for x in 1:nstrains]...)
        tmp = DataFrames.by(tmp, [:variable, :population], :value => x -> round.(x, digits=8))
        unique!(tmp)
        tmp = DataFrames.by(tmp, [:variable, :population], :value_function => x -> x[is_local_minima(x)])
        tmp.chi = [params[1][3] for x in 1:size(tmp, 1)]
        tmp.gamma = vcat([vcat(params[2]...) for x in 1:nstrains]...)
        tmp.sigma = vcat([params[3] for x in 1:nstrains*npopulations]...)
        tmp.R0 = vcat([params[4] for x in 1:nstrains*npopulations]...)
        return(tmp)
    end
end

function explore_parameter_space(;strainstruct = [3 3],
                                 chi           = 0.02,
                                 maxtime       = 10000.0,
                                 gamma         = 0.01:0.33:0.99,
                                 sigma         = 0.05 .* (10 .^ (1:1:4)),
                                 mu            = 0.05,
                                 rnaught       = [2,4,6])
    params = [[γ1, γ2, σ1, σ2, R0] for γ1 in gamma for σ1 in sigma for γ2 in gamma for σ2 in sigma for R0 in rnaught]
    nstrains = trunc(Int, prod(strainstruct))
    npopulations = 2
    results = pmap(x -> clean_eps_timeseries_output!(
                            runMANTIS(strainstructure=strainstruct,
                                      tmax=maxtime,
                                      tstep=(0.95*maxtime):maxtime,
                                      beta=x[5] .* x[3:4], gamma=x[1:2], sigma=x[3:4],
                                      chi=chi, mu=mu)["timeseries"][:,vcat((npopulations * nstrains + 1):(2 * npopulations * nstrains), end)], x, npopulations, nstrains),
                   params)
    return(convert(Matrix, vcat(results...)))
end

function explore_parameter_space_one_pop(;strainstruct = [2 2],
                                          chi           = 0,
                                          maxtime       = 10000.0,
                                          tstep         = 1,
                                          gamma         = 0.02:0.04:0.99,
                                          sigma         = 10.0 .^ (-0.83:0.17:3.25),
                                          mu            = [0.02 0.05 0.1],
                                          rnaught       = [2 4 6])
    params = [[γ, σ, R0, μ] for γ in gamma for σ in sigma for R0 in rnaught for μ in mu]
    nstrains = trunc(Int, prod(strainstruct))
    results = pmap(x -> clean_eps1p_timeseries_output!(
                            runMANTIS(strainstructure=strainstruct,
                                      tmax=maxtime,
                                      tstep=(0.95*maxtime):tstep:maxtime,
                                      beta=x[3] * x[2], gamma=x[1], sigma=x[2],
                                      chi=chi, mu=x[4])["timeseries"][:,vcat((nstrains + 1):(2 * nstrains), end)], x, nstrains),
                   params)
    return(convert(Matrix, vcat(results...)))
end

function variable_movement_rate(;strainstruct = [3 3],
                                 chi           = 10.0 .^ -(1:6),
                                 maxtime       = 10000.0,
                                 gamma         = [0.25 0.75],
                                 sigma         = [1,5,10,50],
                                 mu            = 0.05,
                                 rnaught       = [2,4,6])
    params = [([-χ χ; 0 0], [γ1 γ2], σ, R0) for χ in chi for γ1 in gamma for γ2 in gamma for σ in sigma for R0 in rnaught for rep in 1:10]
    params = params[[x[2][1] != x[2][2] for x in params]]
    nstrains = trunc(Int, prod(strainstruct))
    npopulations = 2
    results = pmap(x -> clean_vmr_timeseries_output!(
                            runMANTIS(strainstructure=strainstruct,
                                      tmax=maxtime,
                                      tstep=(0.95*maxtime):maxtime,
                                      beta=x[4] * x[3], gamma=x[2], sigma=x[3],
                                      chi=x[1], mu=mu)["timeseries"][:,vcat((npopulations * nstrains + 1):(2 * npopulations * nstrains), end)], x, npopulations, nstrains),
                   params)
    return(convert(Matrix, vcat(results...)))
end
