# spectral statistics
function integrate(y::AbstractVector, x::AbstractVector, method::AbstractIntegralAlgorithm=TrapezoidalRule())
    sol = solve(SampledIntegralProblem(y, x), method)
    sol.retcode ≠ ReturnCode.Success && error("solution unsuccessful with code: $(sol.retcode)")
    return sol.u
end

function integrate(y::AbstractMatrix, x::AbstractVector, method::AbstractIntegralAlgorithm=TrapezoidalRule())
    sol = solve(SampledIntegralProblem(y, x; dim=2), method)
    sol.retcode ≠ ReturnCode.Success && error("solution unsuccessful with code: $(sol.retcode)")
    return sol.u
end

function integrate(y::AbstractVector, x::AbstractVector, a::Number, b::Number, method::AbstractIntegralAlgorithm=QuadGKJL(); abstol=nothing)
    f = extrapolate(interpolate((x,), y, Gridded(Linear())), 0.0)
    isnothing(abstol) && (abstol=1.0e-8*unit(eltype(x))*unit(eltype(y)))
    sol = solve(IntegralProblem((x, p)->f(x), a, b), method; abstol)
    sol.retcode ≠ ReturnCode.Success && error("solution unsuccessful with code: $(sol.retcode)")
    return sol.u
end

function integrate(y::AbstractMatrix, x::AbstractVector, a::Number, b::Number, method::AbstractIntegralAlgorithm=QuadGKJL(); abstol=nothing)
    isnothing(abstol) && (abstol=1.0e-8*unit(eltype(x))*unit(eltype(y)))
    sol = ones(eltype(x[1]*y[1,1]), size(y)[1])
    for (i, yᵢ) in enumerate(eachrow(y))
        fᵢ = extrapolate(interpolate((x,), yᵢ, Gridded(Linear())), 0.0)
        solᵢ = solve(IntegralProblem((x, p)->fᵢ(x), a, b), method; abstol)
        solᵢ.retcode ≠ ReturnCode.Success && error("solution $i unsuccessful with code: $(sol.retcode)")
        sol[i] = solᵢ.u
    end
    return sol
end

function integrate(y::AbstractMatrix, x::AbstractVector, a::Union{Number, AbstractVector}, b::Union{Number, AbstractVector}, method::AbstractIntegralAlgorithm=QuadGKJL(); abstol=nothing)
    isnothing(abstol) && (abstol=1.0e-8*unit(eltype(x))*unit(eltype(y)))
    sol = ones(eltype(x[1]*y[1,1]), size(y)[1])
    a = (length(a)≠1) ? a : ones(size(y)[1])*a
    b = (length(b)≠1) ? b : ones(size(y)[1])*b
    for (i, yᵢ) in enumerate(eachrow(y))
        fᵢ = extrapolate(interpolate((x,), yᵢ, Gridded(Linear())), 0.0)
        solᵢ = solve(IntegralProblem((x, p)->fᵢ(x), a[i], b[i]), method; abstol)
        solᵢ.retcode ≠ ReturnCode.Success && error("solution $i unsuccessful with code: $(sol.retcode)")
        sol[i] = solᵢ.u
    end
    return sol
end

∫ = integrate

function spectral_moment(S::AbstractVector, f::AbstractVector, n::Integer; method::AbstractIntegralAlgorithm=TrapezoidalRule())
    return ∫((S.*(f.^n)), f, method)
end

function spectral_moment(S::AbstractMatrix, f::AbstractVector, n::Integer; method::AbstractIntegralAlgorithm=TrapezoidalRule())
    return ∫((S.*(f.^n)'), f, method)
end

function spectral_moment(S::AbstractVector, f::AbstractVector, n::Integer, a::Number, b::Number; method::AbstractIntegralAlgorithm=QuadGKJL())
    return ∫((S.*(f.^n)), f, a, b, method)
end

function spectral_moment(S::AbstractMatrix, f::AbstractVector, n::Integer, a::Union{Number, AbstractVector}, b::Union{Number, AbstractVector}; method::AbstractIntegralAlgorithm=QuadGKJL())
    return ∫((S.*(f.^n)'), f, a, b, method)
end

function energy_period(S::AbstractVecOrMat, f::AbstractVector; method::AbstractIntegralAlgorithm=TrapezoidalRule())
    m_n1 = spectral_moment(S, f, -1; method)
    m_0 = spectral_moment(S, f, 0; method)
    return upreferred.(m_n1./m_0)
end

function energy_period(S::spectrum; method::AbstractIntegralAlgorithm=TrapezoidalRule())
    m_n1 = spectral_moment(S.values, S.frequencies, -1; method)
    m_0 = spectral_moment(S.values, S.frequencies, 0; method)
    return upreferred.(m_n1./m_0)
end

function significant_waveheight(S::AbstractVecOrMat, f::AbstractVector; method::AbstractIntegralAlgorithm=TrapezoidalRule())
	m_0 = spectral_moment(S, f, 0; method)
    return 4(.√m_0)
end

function significant_waveheight(S::spectrum; method::AbstractIntegralAlgorithm=TrapezoidalRule())
    m_0 = spectral_moment(S.values, S.frequencies, 0; method)
    return 4(.√m_0)
end
