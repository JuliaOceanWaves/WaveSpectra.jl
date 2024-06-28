# WaveSpectra.jl Documentation

```@contents
```

```@meta
DocTestSetup = quote
    function box(x)
        ax = ustrip(x)
        if (ax == 0.5 || ax == 1.5)
            return 0.5m^2
        elseif (0.5 < ax < 1.5)
            return 1.0m^2
        else
            return 0.0m^2
        end
    end
    f = range(0.1, 2.1, 51)
end
end
```

## Functions
```@docs
DiscreteOmnidirectionalSpectrum
OmnidirectionalSpectrum
frequency_unit
frequency_dimension
quantity
spectral_moment
energy_period
significant_waveheight
convert_frequency
```

## Index

```@index
```
