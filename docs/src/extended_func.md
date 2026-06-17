# Extended Functionality

## Dispersion Relations


The DispersionRelation.jl submodule contains the dispersion relations for linear interfacial
waves. Currently only gravity wave dispersion relations are implemented.

!!! note

    We use phase velocity ($V_p$)

See also [https://en.wikipedia.org/wiki/Dispersion\_(water_waves)](https://en.wikipedia.org/wiki/Dispersion_(water_waves)) 
for more information on dispersion.

```julia

julia> x = (1.0:3:10.0) * Hz;

julia> S = OmnidirectionalSpectrum(([1.0, 6.0, 3.0, 2.0])*m^2/Hz, x)
4-element OmnidirectionalSpectrum{m² Hz⁻¹}{Hz}
Spectral density of the quantity (m²):
  • Axis: Frequency (Hz)
and data(m² Hz⁻¹):
 1.0
 6.0
 3.0
 2.0

julia> uconvert(rad/m, :axis, S, WaveSpectra.DispersionRelations.gravitywaves_deepwater())
4-element OmnidirectionalSpectrum{m³ rad⁻¹}{rad m⁻¹}
Spectral density of the quantity (m²):
  • Axis: Angular Wavenumber (rad m⁻¹)
and data(m³ rad⁻¹):
 0.12420267319576646
 0.1863040097936497
 0.053229717083899904
 0.02484053463915329

julia> uconvert(rad/m, :axis, S, WaveSpectra.DispersionRelations.gravitywaves_shallowwater(3m))
4-element OmnidirectionalSpectrum{m³ rad⁻¹}{rad m⁻¹}
Spectral density of the quantity (m²):
  • Axis: Angular Wavenumber (rad m⁻¹)
and data(m³ rad⁻¹):
 0.863258964143784
 5.179553784862704
 2.589776892431352
 1.726517928287568

julia> uconvert(rad/m, :axis, S, WaveSpectra.DispersionRelations.gravitywaves(3m))
4-element OmnidirectionalSpectrum{m³ rad⁻¹}{rad m⁻¹}
Spectral density of the quantity (m²):
  • Axis: Angular Wavenumber (rad m⁻¹)
and data(m³ rad⁻¹):
 0.12420267338189336
 0.1863040097936497
 0.053229717083899904
 0.02484053463915329

```

Please refer to the full syntax for each function [here](@ref dispersion_relation_syntax).


## [Parametric Spectra](@id parametric_spectra) 


Brief description of parametric spectrum

```julia
julia> f = (1.0:0.5:10.0) .* Hz;

julia> PM = WaveSpectra.ParametricSpectra.spectrum_pierson_moskowitz(f, 8.0 * m, 4.0 * Hz, peak_frequency=true)
19-element OmnidirectionalSpectrum{m² Hz⁻¹}{Hz}
Spectral density of the quantity (m²):
  • Axis: Frequency (Hz)
and data(m² Hz⁻¹):
 1.3423912370794646e-72
 4.849071845233627e-13
 0.001701611389553719
 0.335309891337633
 1.3421276585162367
 1.6633580565889663
 ⋮
 0.08116740160851418
 0.0604917892588704
 0.045764352079589225
 0.035103816684478595
 0.02727015323859412

```

```julia
julia> f = (1.0:0.5:10.0) .* Hz;

julia> JS = WaveSpectra.ParametricSpectra.spectrum_jonswap(f, 8.0 * m, 4.0 * Hz, 3)
19-element OmnidirectionalSpectrum{m² Hz⁻¹}{Hz}
Spectral density of the quantity (m²):
  • Axis: Frequency (Hz)
and data(m² Hz⁻¹):
 1.9250099830902135e-88
 3.0924699228600333e-16
 0.00014581161763694455
 0.10965073255059529
 0.757687377655647
 2.989952938117487
 1.6476870418863625
 0.7639004421371922
 0.5256674689081673
 0.3629932785530433
 0.2513107958135555
 0.17605740404750989
 0.1252620883911322
 0.09060488782908652
 0.06661272946290435
 0.04973940996324117
 0.03768327349740816
 0.02893628972068863
 0.022497539992942253
```

```julia
julia> Θ = (0:15:90)*°; f = (1:1:10)*Hz;

julia> spread_func = WaveSpectra.ParametricSpectra.spread_cartwright(Θ, f, 10*°, 20*°);

julia> S_omni = OmnidirectionalSpectrum((1:10)*m, f);

julia> S = Spectrum(S_omni, spread_func)
10×7 Spectrum{m °⁻¹}{Hz}{°}
Spectral density of the quantity (Hz m) with polar coordinates:
  • Axis 1: Frequency (Hz)
  • Axis 2: Direction (°)
and data(m °⁻¹):
 0.025494124746630847  0.027844248353181597  0.01788639292257555   0.006653152827856875  0.0013815841198168257  0.00015035145677660058  7.749226287575886e-6
 0.050988249493261695  0.055688496706363194  0.0357727858451511    0.01330630565571375   0.0027631682396336513  0.00030070291355320115  1.5498452575151773e-5
 0.07648237423989254   0.08353274505954479   0.053659178767726655  0.019959458483570624  0.004144752359450477   0.00045105437032980176  2.324767886272766e-5
 0.10197649898652339   0.11137699341272639   0.0715455716903022    0.0266126113114275    0.005526336479267303   0.0006014058271064023   3.0996905150303545e-5
 0.12747062373315424   0.139221241765908     0.08943196461287775   0.03326576413928437   0.006907920599084128   0.0007517572838830029   3.874613143787943e-5
 0.15296474847978508   0.16706549011908958   0.10731835753545331   0.03991891696714125   0.008289504718900954   0.0009021087406596035   4.649535772545532e-5
 0.17845887322641593   0.19490973847227117   0.12520475045802887   0.046572069794998124  0.00967108883871778    0.001052460197436204    5.4244584013031204e-5
 0.20395299797304678   0.22275398682545278   0.1430911433806044    0.053225222622855     0.011052672958534605   0.0012028116542128046   6.199381030060709e-5
 0.22944712271967763   0.2505982351786344    0.16097753630317996   0.059878375450711875  0.01243425707835143    0.0013531631109894053   6.974303658818298e-5
 0.2549412474663085    0.278442483531816     0.1788639292257555    0.06653152827856874   0.013815841198168257   0.0015035145677660057   7.749226287575886e-5
```

Please refer to the full syntax for each function [here](@ref parametric_spectra_syntax).


## [Moments](@id moments)


The following examples are different functions used in literature for characterizing 
spectra. 

```julia
julia> x = (1.0:3:10.0) * Hz; S = OmnidirectionalSpectrum(([1.0, 6.0, 3.0, 2.0])*m^2, x);

julia> WaveSpectra.Moments.moment(S, 0)
31.5 m²

julia> WaveSpectra.Moments.energy_frequency(S)
4.152542372881356 Hz

julia> WaveSpectra.Moments.mean_frequency(S)
5.285714285714286 Hz

julia> WaveSpectra.Moments.mean_wavelength(S)
0.09051335476420994 m

julia> WaveSpectra.Moments.significant_waveheight(S)
22.44994432064365 m

julia> WaveSpectra.Moments.steepness(S)
248.02908232852977

julia> WaveSpectra.Moments.zero_crossing_frequency(S)
5.719640348333601 Hz

```

Please refer to the full syntax for each function [here](@ref moments_syntax).

## [Spectral Shapes](@id spectral_shapes)

Spectral Shape

```julia
julia> S = [1.0, 2.0, 3.5, 1.0, 0.5]*m^2; f = (1.0:1.0:5.0)*Hz;

julia> S_omni = OmnidirectionalSpectrum(S, f);

julia> S_shape = WaveSpectra.Shapes.OmnidirectionalSpectrumShape(S_omni)
5-element OmnidirectionalSpectrum{1}{}
Spectral density of the quantity ():
  • Axis: Ndfrequency ()
and data(1):
 0.02106741573033708
 0.04213483146067416
 0.07373595505617979
 0.02106741573033708
 0.01053370786516854

 julia> S_omni_scaled = OmnidirectionalSpectrum(S_shape, 5m, 2Hz)
5-element OmnidirectionalSpectrum{m² Hz⁻¹}{Hz}
Spectral density of the quantity (m²):
  • Axis: Frequency (Hz)
and data(m² Hz⁻¹):
 0.2633426966292135
 0.526685393258427
 0.9216994382022473
 0.2633426966292135
 0.13167134831460675

julia> S_omni_scaled = WaveSpectra.Shapes.scale(S_omni, 5m, 2Hz)
5-element OmnidirectionalSpectrum{m² Hz⁻¹}{Hz}
Spectral density of the quantity (m²):
  • Axis: Frequency (Hz)
and data(m² Hz⁻¹):
 0.2633426966292135
 0.526685393258427
 0.9216994382022473
 0.2633426966292135
 0.13167134831460675
```


Please refer to the full syntax for each function [here](@ref spectral_shapes_syntax).


## Syntax

  - [Dispersion Relation](@ref dispersion_relation_syntax)
  - [Parametric Spectra](@ref parametric_spectra_syntax)
  - [Moments](@ref moments_syntax)
  - [Spectral Shapes](@ref spectral_shapes_syntax)


### [Dispersion Relation](@id dispersion_relation_syntax)

```@autodocs; canonical=false
Modules = [WaveSpectra.DispersionRelations]
Order = [:function]
```

### [Parametric Spectra](@id parametric_spectra_syntax)

```@autodocs; canonical=false
Modules = [WaveSpectra.ParametricSpectra]
Order = [:function, :type]
```

### [Moments](@id moments_syntax)

```@autodocs; canonical=false
Modules = [WaveSpectra.Moments]
Order = [:function]
```

### [Spectral Shape](@id spectral_shapes_syntax)

```@autodocs; canonical=false
Modules = [WaveSpectra.Shapes]
Order = [:function, :type]
```
