# [Other Capabilities](@id other_funcs)

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

## [Spectral Shape](@id spectral_shape)

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


Please refer to the full syntax for each function [here](@ref spectral_shape_syntax).

## Syntax

### [Parametric Spectra](@id parametric_spectra_syntax)

```@autodocs; canonical=false
Modules = [WaveSpectra.ParametricSpectra]
Order = [:function, :type]
```

### [Spectral Shape](@id spectral_shape_syntax)

```@autodocs; canonical=false
Modules = [WaveSpectra.Shapes]
Order = [:function, :type]
```