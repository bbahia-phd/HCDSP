export SeisAddNoise,
        MeasureSNR,
        Hamming,
        add_erratic_noise,
        decimate_traces,
        linear_dispersive_events,
        irregular_hyperbolic_events,
        SeisLinearEvents,
        SeisHypEvents,
        SeisParabEvents


    include("add_random_noise.jl")
    include("add_erratic_noise.jl")
    include("decimate_traces.jl")
    include("linear_dispersive_events.jl")
    include("SeisHypEvents.jl")
    include("SeisParabEvents.jl")
    include("SeisLinearEvents.jl")