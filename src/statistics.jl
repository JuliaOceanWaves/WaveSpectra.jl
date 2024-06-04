

function integrate(func::Function, f_begin::Union{Number, Nothing}=nothing, f_end::Union{Number, Nothing}=nothing; alg::AbstractIntegralAlgorithm=QuadGKJL())
    isnothing(f_begin) && (f_begin = isunitful(func) ? 0Hz : 0)
    isnothing(f_end) && (f_end = isunitful(func) ? Inf*Hz : Inf)
    sol = solve(IntegralProblem((x,p)->func(x), (f_begin, f_end), nothing), alg)
    sol.retcode ≠ ReturnCode.Success && error("solution unsuccessful with code: $(sol.retcode)")
    return sol.u
end

function integrate(S::WaveSpectrum, f_begin::Union{Number, Nothing}=nothing, f_end::Union{Number, Nothing}=nothing; alg::AbstractIntegralAlgorithm=QuadGKJL())
    isnothing(f_begin) && (f_begin = isunitful(S) ? 0Hz : 0)
    isnothing(f_end) && (f_end = isunitful(S) ? Inf*Hz : Inf)
    return integrate(S.func, f_begin, f_end; alg=alg)
end

∫ = integrate

function spectral_moment(S::WaveSpectrum, n::Int, f_begin::Union{Number, Nothing}=nothing, f_end::Union{Number, Nothing}=nothing; alg::AbstractIntegralAlgorithm=QuadGKJL())
    return ∫(f -> S(f)*f^n, f_begin, f_end; alg=alg)
end

function energy_period(S::WaveSpectrum, f_begin::Union{Number, Nothing}=nothing, f_end::Union{Number, Nothing}=nothing; alg::AbstractIntegralAlgorithm=QuadGKJL())
    m_n1 = spectral_moment(S, -1, f_begin, f_end; alg=alg)
    m_0 = spectral_moment(S, 0, f_begin, f_end; alg=alg)
    return upreferred(m_n1/m_0)
end

function significant_waveheight(S::WaveSpectrum, f_begin::Union{Number, Nothing}=nothing, f_end::Union{Number, Nothing}=nothing; alg::AbstractIntegralAlgorithm=QuadGKJL())
	m_0 = spectral_moment(S, 0, f_begin, f_end; alg=alg)
    return 4√m_0
end
