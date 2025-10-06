# cosmology.jl

A Julia package for calculating cosmological quantities including distances, Hubble parameters, and linear growth factors.

## Overview

`cosmology.jl` provides tools for computing standard cosmological quantities for ΛCDM and wCDM cosmologies, with support for massive neutrinos and radiation. Originally implemented in Python (2013), translated to Julia and packaged (2017-2019).

## Features

- **Cosmological distances**: comoving, angular diameter, luminosity distances
- **Hubble parameter** H(z) with full matter/radiation/dark energy evolution
- **Linear growth factors** D(z) and dlnD/dlna via numerical integration
- **Massive neutrino support** with arbitrary neutrino masses
- **Redshift-distance conversions** via spline interpolation
- **Volume calculations** for survey geometry
- **Pre-computed tables** cached for performance

## Installation

### From GitHub
```julia
using Pkg
Pkg.add(url="https://github.com/yourusername/cosmology.jl")
```

### Development Installation
```julia
Pkg.develop(path="/path/to/cosmology")
```

## Quick Start

```julia
using cosmology

# Initialize cosmology (Planck 2018-like parameters)
cosmo = cosParams(
    h = 0.6766,        # Hubble constant / 100
    Om = 0.3111,       # Matter density parameter
    Ode = 0.6889,      # Dark energy density parameter
    Ob = 0.048975,     # Baryon density parameter
    wde = -1.0,        # Dark energy equation of state
    unit = "Mpc/h"     # Distance units: "Mpc/h" or "Mpc"
)

# Calculate cosmological quantities
z = 2.7
H_z = Hz(cosmo, z)              # Hubble parameter at redshift z
chi = chiz(cosmo, z)            # Comoving distance
dA = dAz(cosmo, z)              # Angular diameter distance
dL = dLz(cosmo, z)              # Luminosity distance
D_growth = Dplusz(cosmo, z)     # Linear growth factor D(z)/D(0)
f_growth = dlnDdlna(cosmo, z)   # dlnD/dlna (growth rate)
```

## API Reference

### Initialization

#### `cosParams(h, Om, Ode; kwargs...)`

Create a cosmological parameter object.

**Required arguments:**
- `h`: Hubble constant in units of 100 km/s/Mpc
- `Om`: Matter density parameter at z=0
- `Ode`: Dark energy density parameter at z=0

**Optional keyword arguments:**
- `Ob=0.048555`: Baryon density parameter
- `wde=-1.0`: Dark energy equation of state parameter
- `unit="Mpc/h"`: Distance units ("Mpc/h" or "Mpc")
- `Tcmb=2.726`: CMB temperature (K), set to 0 to ignore radiation
- `mnu=[0.06, 0., 0.]`: Neutrino masses in eV (3 species)
- `calc_fxnu=false`: Force recalculation of neutrino energy density
- `calc_growth=false`: Force recalculation of growth factors
- `calc_chiz=false`: Force recalculation of distance-redshift relation

**Returns:** `cosParams` struct with all cosmological functions

### Distance Functions

- `chiz(cosmo, z)`: Comoving distance (line-of-sight)
- `dAz(cosmo, z)`: Angular diameter distance
- `dAcz(cosmo, z)`: Comoving angular diameter distance = (1+z)·dA(z)
- `dLz(cosmo, z)`: Luminosity distance = (1+z)²·dA(z)
- `findz(cosmo, chi)`: Find redshift corresponding to comoving distance chi

### Hubble Parameter

- `Hz(cosmo, z)`: Hubble parameter H(z)
- `Ha(cosmo, a)`: Hubble parameter H(a) as function of scale factor
- `aHa(cosmo, a)`: a·H(a)
- `asqHa(cosmo, a)`: a²·H(a)

### Density Parameters

- `Omz(cosmo, z)`: Matter density parameter Ωₘ(z)
- `Odez(cosmo, z)`: Dark energy density parameter ΩΛ(z)
- `Orz(cosmo, z)`: Radiation density parameter Ωᵣ(z)
- `Okz(cosmo, z)`: Curvature density parameter Ωₖ(z)

### Growth Functions

- `Dz(cosmo, z)`: Linear growth factor D(z)
- `Dplusz(cosmo, z)`: Normalized growth D(z)/D(0)
- `dlnDdlna(cosmo, z)`: Logarithmic derivative dlnD/dlna

### Other Functions

- `agez(cosmo, z)`: Age of universe at redshift z
- `Volume(cosmo, zi, zf, Area)`: Comoving volume between zi and zf for solid angle Area (deg²)
- `comPtlHorizon(cosmo, z)`: Comoving particle horizon

## Examples

### Example 1: Standard ΛCDM Cosmology (Planck 2018)

```julia
using cosmology

cosmo = cosParams(
    h = 0.6766,
    Om = 0.3111,
    Ode = 0.6889,
    Ob = 0.048975
)

# Distance modulus for supernova cosmology
z_sn = 0.5
dL = dLz(cosmo, z_sn)
μ = 5 * log10(dL) + 25  # Distance modulus (if dL in Mpc)
```

### Example 2: Growth Rate for RSD

```julia
# Calculate f·σ₈ for redshift space distortions
z = 1.0
f = dlnDdlna(cosmo, z) - 1  # Linear growth rate
D_ratio = Dplusz(cosmo, z)
# f·σ₈(z) ≈ f(z) · σ₈(0) · D(z)/D(0)
```

### Example 3: Survey Volume

```julia
# Calculate survey volume for HETDEX
z_min, z_max = 1.9, 3.5
survey_area = 420.0  # deg²

V = Volume(cosmo, z_min, z_max, survey_area)
println("Survey volume: $V (Mpc/h)³")
```

### Example 4: Massive Neutrinos

```julia
# Cosmology with 0.06 eV neutrino mass
cosmo_nu = cosParams(
    h = 0.6766,
    Om = 0.3111,
    Ode = 0.6889,
    mnu = [0.06, 0.0, 0.0]  # One massive neutrino
)

# Neutrino mass affects H(z) and distances at high-z
println("Ωₘ(z=0) = $(Omz(cosmo_nu, 0.0))")
println("Ωₘ(z=3) = $(Omz(cosmo_nu, 3.0))")
```

## Performance Notes

- **First run**: Calculates and caches lookup tables in `./tables/` directory
  - Neutrino energy density function `frhonu_xnu.dat`
  - Growth factors `gfactor_Om=*_Ode=*_wde=*.dat`
  - Distance-redshift relation `chiz_*.dat`

- **Subsequent runs**: Loads pre-computed tables instantly via spline interpolation

- To force recalculation, use `calc_fxnu=true`, `calc_growth=true`, `calc_chiz=true` or delete the `./tables/` directory

## Physical Constants

The package uses the following constants:
- Hubble constant: H₀ = 100h km/s/Mpc = 0.00033356h Mpc⁻¹
- Critical density: ρcrit = 2.7751×10¹¹ h² M⊙ (Mpc/h)⁻³
- Speed of light: c = 299792.458 km/s
- CMB temperature (pivot): Tcmb = 2.726 K

## Citation

If you use this package in your research, please cite:

```bibtex
@software{cosmology_jl,
  author = {Donghui Jeong},
  title = {cosmology.jl: Cosmological Calculations in Julia},
  year = {2019},
  url = {https://github.com/yourusername/cosmology.jl}
}
```

## References

For equations and methodology, see `cosmo_eqns.pdf` (if available in repository).

## Requirements

- Julia ≥ 1.6
- Dierckx.jl (spline interpolation)
- OrdinaryDiffEq.jl (growth factor integration)
- QuadGK.jl (numerical integration)
- SpecialFunctions.jl (special functions)

## License

See LICENSE file for details.

## Author

Donghui Jeong
Pennsylvania State University
djeong@psu.edu

---

**Original Python implementation**: March 2013
**Julia translation**: October 2017
**Package release**: January 2019
