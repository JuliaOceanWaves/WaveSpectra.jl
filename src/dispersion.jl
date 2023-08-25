# TODO: Add all other waves in https://en.wikipedia.org/wiki/Dispersion_(water_waves).
struct DeepWater{g, } <: Equivalence end

edconvert(::dimtype(AngularVelocity), k::AngularWavenumber, ::DeepWater) = √(g*k)
edconvert(::dimtype(AngularWavenumber), ω::AngularVelocity, ::DeepWater) = ω^2/g


function expand_dispersion(d::Equivalence)
    # d defines ω <=> k and nothing else
    function edconvert()
