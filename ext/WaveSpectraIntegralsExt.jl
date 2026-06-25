module WaveSpectraIntegralsExt

using Integrals
using WaveSpectra

function Integrals.solve(
        prob::Integrals.SampledIntegralProblem,
        alg::WaveSpectra.RectangularRule;
        kwargs...
)
    NamedTuple(kwargs) == NamedTuple() ||
        throw(ArgumentError("There are no keyword arguments allowed to `solve`"))

    u = prob.y isa AbstractVector ?
        WaveSpectra._solve_sampled_integral(prob.y, prob.x; method = alg) :
        WaveSpectra._solve_sampled_integral(prob.y, prob.x; method = alg, dim = prob.dim)
    return (u = u,)
end

end
