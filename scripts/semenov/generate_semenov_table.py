#!/usr/bin/env python3
"""
generate_semenov_table.py

Generate a dense Rosseland-mean opacity table for the Semenov et al. (2003)
"normal silicate / homogeneous sphere" model (nrm / h / s), using only Python.

The script reproduces the logic of the original FORTRAN 77 package:

  - dust-dominated opacity:
      fifth-degree polynomial fits stored in subroutine nrm_h_s()
  - density-dependent evaporation temperatures:
      logic from subroutine cop()
  - smoothing across evaporation transitions:
      logic from subroutine cop()
  - gas-dominated opacity:
      kR_h2001.dat + interpolation logic from gop()/EINT()

No FORTRAN compiler or runtime is required.

Inputs
------
1. opacity.f
2. kR_h2001.dat
3. optionally kR.out, used only as a regression check against the original
   FORTRAN executable output.

Outputs
-------
1. semenov_rosseland_nrm_h_s.npz
   Binary NumPy file containing:
       logrho      [NRHO]
       logT        [NTEMP]
       logkappa    [NRHO, NTEMP]
       rho_cgs     [NRHO]
       T_K         [NTEMP]
       kappa_cgs   [NRHO, NTEMP]

2. semenov_rosseland_nrm_h_s.dat
   Plain ASCII table suitable for later conversion/reading by FARGO3D.

3. fig_semenov_rosseland_table.png
   Diagnostic opacity map.

Notes
-----
- rho is in g cm^-3.
- T is in K.
- kappa_R is in cm^2 g^-1.
- The gas-opacity data in kR_h2001.dat are log10(kappa_R).
- The original EINT() interpolation is linear in T and rho for the tabulated
  log10(kappa), with a quadratic low-density extrapolation permitted by gop().
- The dense production table is sampled uniformly in log10(rho) and log10(T).
  Later FARGO3D interpolation should be bilinear in (logrho, logT, logkappa).
"""

import argparse
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# SEMENOV MODEL CONSTANTS FROM cop()
# ============================================================

# Density knots used for density-dependent evaporation temperatures.
EVAP_RHO = np.array(
    [1.0e-18, 1.0e-16, 1.0e-14, 1.0e-12,
     1.0e-10, 1.0e-08, 1.0e-06, 1.0e-04],
    dtype=np.float64,
)

# Columns:
#   0 water ice
#   1 metallic iron
#   2 orthopyroxene
#   3 olivine
EVAP_T = np.array(
    [
        [109.0,  835.0,  902.0,  929.0],
        [118.0,  908.0,  980.0,  997.0],
        [129.0,  994.0, 1049.0, 1076.0],
        [143.0, 1100.0, 1129.0, 1168.0],
        [159.0, 1230.0, 1222.0, 1277.0],
        [180.0, 1395.0, 1331.0, 1408.0],
        [207.0, 1612.0, 1462.0, 1570.0],
        [244.0, 1908.0, 1621.0, 1774.0],
    ],
    dtype=np.float64,
)

# Widths used by COP() to smooth opacity changes around evaporation boundaries.
SMOOTH_DT = np.array([5.0, 5.0, 15.0, 5.0, 100.0], dtype=np.float64)

# Hard validity checks used by the original gas-opacity routine gop().
GAS_RHO_MIN = 1.0e-19
GAS_RHO_MAX = 1.0e-07
GAS_T_MIN = 500.0
GAS_T_MAX = 10000.0


# ============================================================
# FORTRAN SOURCE PARSING
# ============================================================

_NUM_RE = re.compile(
    r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[DEde][-+]?\d+)?"
)


def _numbers(text):
    """Extract FORTRAN-style real numbers from text."""
    values = []
    for line in text.splitlines():
        # FORTRAN source uses ! for inline comments in these DATA blocks.
        line = line.split("!", 1)[0]
        for token in _NUM_RE.findall(line):
            values.append(float(token.replace("D", "E").replace("d", "e")))
    return np.asarray(values, dtype=np.float64)


def parse_nrm_h_s_rosseland_coefficients(opacity_f):
    """
    Read the Rosseland dust-fit coefficients directly from nrm_h_s().

    FORTRAN stores each region as:
        eR(i,6), eR(i,5), ..., eR(i,1)

    while COP evaluates:
        ((((eR(i,1) T + eR(i,2)) T + ...) T + eR(i,6))

    Therefore each six-value DATA block must be reversed before np.polyval().
    """
    source = Path(opacity_f).read_text(encoding="latin-1")

    section = re.search(
        r"SUBROUTINE\s+nrm_h_s\b(.*?)SUBROUTINE\s+nrm_h_a\b",
        source,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if section is None:
        raise RuntimeError("Could not locate SUBROUTINE nrm_h_s in opacity.f.")

    rosseland = re.search(
        r"DATA\s+eR\(1,6\).*?&\s*/(.*?)/",
        section.group(1),
        flags=re.IGNORECASE | re.DOTALL,
    )
    if rosseland is None:
        raise RuntimeError("Could not locate the nrm_h_s Rosseland DATA block.")

    values = _numbers(rosseland.group(1))

    if values.size != 30:
        raise RuntimeError(
            f"Expected 30 nrm_h_s Rosseland coefficients, found {values.size}."
        )

    # Text order per region: constant, T, T^2, ..., T^5.
    text_order = values.reshape(5, 6)

    # np.polyval order: T^5, T^4, ..., constant.
    return text_order[:, ::-1].copy()


def parse_gop_grids(opacity_f):
    """Extract the original 71-point gas temperature and density grids."""
    source = Path(opacity_f).read_text(encoding="latin-1")

    section = re.search(
        r"SUBROUTINE\s+gop\b(.*?)SUBROUTINE\s+bint\b",
        source,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if section is None:
        raise RuntimeError("Could not locate SUBROUTINE gop in opacity.f.")

    body = section.group(1)

    tmatch = re.search(
        r"DATA\s+T\(1\).*?&\s*/(.*?)/",
        body,
        flags=re.IGNORECASE | re.DOTALL,
    )
    rmatch = re.search(
        r"DATA\s+rho\(1\).*?&\s*/(.*?)/",
        body,
        flags=re.IGNORECASE | re.DOTALL,
    )

    if tmatch is None or rmatch is None:
        raise RuntimeError("Could not parse the gas T/rho grids from gop().")

    T_grid = _numbers(tmatch.group(1))
    rho_grid = _numbers(rmatch.group(1))

    if T_grid.size != 71 or rho_grid.size != 71:
        raise RuntimeError(
            "Expected 71 gas temperature and density knots; "
            f"found NT={T_grid.size}, NRHO={rho_grid.size}."
        )

    if not np.all(np.diff(T_grid) > 0.0):
        raise RuntimeError("Semenov gas temperature grid is not increasing.")

    if not np.all(np.diff(rho_grid) < 0.0):
        raise RuntimeError("Semenov gas density grid is not decreasing.")

    return rho_grid, T_grid


def read_kR_h2001(path):
    """
    Read kR_h2001.dat exactly as init_g() does.

    The first scalar is a dummy value. The following 71*71 values are
    log10(kappa_R), and FORTRAN fills eG(i,j) with i (rho) outermost and
    j (temperature) innermost.
    """
    text = Path(path).read_text(encoding="latin-1")
    values = _numbers(text)

    if values.size != 1 + 71 * 71:
        raise RuntimeError(
            f"Expected 5042 numbers in kR_h2001.dat, found {values.size}."
        )

    dummy = values[0]
    eG = values[1:].reshape(71, 71)

    return dummy, eG


# ============================================================
# ORIGINAL SEMENOV INTERPOLATION LOGIC
# ============================================================

def bint_linear_extrap(xa, x, y):
    """
    Reproduce bint() from opacity.f.

    Despite the FORTRAN comment calling this 'quadratic interpolation',
    the implementation is piecewise linear interpolation with linear
    extrapolation using the first or last interval.
    """
    x = np.asarray(x)
    y = np.asarray(y)

    if xa < x[0]:
        i2 = 1
    elif xa > x[-1]:
        i2 = x.size - 1
    else:
        i2 = int(np.searchsorted(x, xa, side="left"))
        if i2 == 0:
            i2 = 1

    i1 = i2 - 1

    return y[i1] + (
        (y[i2] - y[i1])
        / (x[i2] - x[i1])
        * (xa - x[i1])
    )


def eint_logkappa(rho_in, T_in, rho_grid, T_grid, eG):
    """
    Reproduce EINT() for the Semenov gas table.

    Interpolation order follows the original FORTRAN:
      1. linear interpolation in T at every rho knot;
      2. linear interpolation in rho;
      3. for rho below the last tabulated knot, quadratic extrapolation
         using the three lowest-density points.

    eG stores log10(kappa_R).
    """
    # Original EINT requires a lower temperature knot and an upper knot.
    it_hi = int(np.searchsorted(T_grid, T_in, side="left"))

    if it_hi <= 0:
        raise ValueError(
            f"EINT cannot interpolate T={T_in:g} K at/below "
            f"its first gas knot {T_grid[0]:g} K."
        )

    if it_hi >= T_grid.size:
        # At exactly the last knot searchsorted returns its index, so this
        # branch is only for truly larger values.
        if np.isclose(T_in, T_grid[-1]):
            it_hi = T_grid.size - 1
        else:
            raise ValueError(
                f"EINT temperature {T_in:g} K exceeds gas grid maximum."
            )

    it_lo = it_hi - 1

    frac_T = (
        (T_in - T_grid[it_lo])
        / (T_grid[it_hi] - T_grid[it_lo])
    )

    xtmp = (
        eG[:, it_lo]
        + frac_T * (eG[:, it_hi] - eG[:, it_lo])
    )

    # rho_grid is monotonically decreasing.
    candidates = np.flatnonzero(rho_in >= rho_grid)

    if candidates.size:
        i_hi = int(candidates[0])
        i_lo = i_hi - 1

        if i_lo >= 0:
            return xtmp[i_lo] + (
                (xtmp[i_hi] - xtmp[i_lo])
                / (rho_grid[i_hi] - rho_grid[i_lo])
                * (rho_in - rho_grid[i_lo])
            )

    # Same quadratic low-density extrapolation as EINT().
    n = rho_grid.size

    a = (
        (xtmp[n-1] - xtmp[n-2])
        / (rho_grid[n-1] - rho_grid[n-2])
        -
        (xtmp[n-2] - xtmp[n-3])
        / (rho_grid[n-2] - rho_grid[n-3])
    ) / (rho_grid[n-1] - rho_grid[n-3])

    b = (
        (xtmp[n-1] - xtmp[n-2])
        / (rho_grid[n-1] - rho_grid[n-2])
        - a * (rho_grid[n-1] + rho_grid[n-2])
    )

    c = (
        xtmp[n-1]
        - a * rho_grid[n-1]**2
        - b * rho_grid[n-1]
    )

    return a * rho_in**2 + b * rho_in + c


def gas_rosseland(
    rho_cgs,
    T_K,
    rho_grid,
    T_grid,
    eG,
):
    """Python reproduction of gop() for Rosseland mean opacity."""
    if rho_cgs > GAS_RHO_MAX or rho_cgs < GAS_RHO_MIN:
        return 0.0

    if T_K < GAS_T_MIN or T_K > GAS_T_MAX:
        return 0.0

    logkappa = eint_logkappa(
        rho_cgs,
        T_K,
        rho_grid,
        T_grid,
        eG,
    )

    return 10.0**logkappa


# ============================================================
# FULL SEMENOV nrm/h/s ROSSeland OPACITY
# ============================================================

def evaporation_boundaries(rho_cgs):
    """
    Reproduce the five transition temperatures assembled inside COP().
    """
    T_ev = np.array(
        [
            bint_linear_extrap(
                rho_cgs,
                EVAP_RHO,
                EVAP_T[:, j],
            )
            for j in range(4)
        ],
        dtype=np.float64,
    )

    boundaries = np.empty(5, dtype=np.float64)

    # Water ice.
    boundaries[0] = T_ev[0]

    # Volatile organics.
    boundaries[1] = 275.0

    # Refractory organics.
    boundaries[2] = 425.0

    # Troilite.
    boundaries[3] = 680.0

    # Iron / orthopyroxene / olivine effective transition.
    tmax1 = max(T_ev[1], T_ev[2])
    tmax2 = max(T_ev[2], T_ev[3])
    boundaries[4] = min(tmax1, tmax2)

    return boundaries


def dust_fit(coefficients, region, T_K):
    """Evaluate one fifth-degree nrm/h/s Rosseland dust fit."""
    return float(np.polyval(coefficients[region], T_K))


def semenov_rosseland(
    rho_cgs,
    T_K,
    dust_coeff,
    gas_rho_grid,
    gas_T_grid,
    eG,
):
    """
    Full Python reproduction of COP() for nrm / h / s Rosseland opacity.

    Returns kappa_R in cm^2 g^-1.
    """
    if T_K < 1.0:
        return 0.0

    boundaries = evaporation_boundaries(rho_cgs)

    # FORTRAN default: KK=6 -> gas.
    region = 5

    if T_K <= boundaries[0] + SMOOTH_DT[0]:
        region = 0

    for it in range(1, 5):
        if (
            T_K > boundaries[it-1] + SMOOTH_DT[it-1]
            and T_K <= boundaries[it] + SMOOTH_DT[it]
        ):
            region = it

    # Gas-dominated region.
    if region == 5:
        return gas_rosseland(
            rho_cgs,
            T_K,
            gas_rho_grid,
            gas_T_grid,
            eG,
        )

    # COP() switches to smoothing if T lies near ANY transition.
    smooth = np.any(
        np.abs(T_K - boundaries) <= SMOOTH_DT
    )

    if not smooth:
        return dust_fit(dust_coeff, region, T_K)

    T1 = boundaries[region] - SMOOTH_DT[region]
    T2 = boundaries[region] + SMOOTH_DT[region]
    TD = T_K - boundaries[region]

    k_left = dust_fit(dust_coeff, region, T1)

    if region == 4:
        k_right = gas_rosseland(
            rho_cgs,
            T2,
            gas_rho_grid,
            gas_T_grid,
            eG,
        )
    else:
        k_right = dust_fit(
            dust_coeff,
            region + 1,
            T2,
        )

    AA = 0.5 * (k_left - k_right)
    BB = 0.5 * (k_left + k_right)
    FF = np.pi / (2.0 * SMOOTH_DT[region])

    return BB - AA * np.sin(FF * TD)


# ============================================================
# REGRESSION AGAINST PROVIDED kR.out
# ============================================================

def validate_against_kR_out(
    kR_out,
    rho_reference,
    dust_coeff,
    gas_rho_grid,
    gas_T_grid,
    eG,
):
    """
    Compare Python COP() reproduction with the provided original FORTRAN
    output file. kR.out itself is printed with limited decimal precision,
    so ~1e-5 relative differences from text rounding are expected.
    """
    path = Path(kR_out)

    if not path.exists():
        return None

    ref = np.loadtxt(path)

    if ref.ndim != 2 or ref.shape[1] < 2:
        raise RuntimeError(f"Unexpected format in {path}.")

    T_ref = ref[:, 0]
    k_ref = ref[:, 1]

    k_py = np.array(
        [
            semenov_rosseland(
                rho_reference,
                T,
                dust_coeff,
                gas_rho_grid,
                gas_T_grid,
                eG,
            )
            for T in T_ref
        ],
        dtype=np.float64,
    )

    valid = k_ref > 0.0

    rel = np.full_like(k_ref, np.nan)
    rel[valid] = (k_py[valid] - k_ref[valid]) / k_ref[valid]

    return {
        "T": T_ref,
        "k_ref": k_ref,
        "k_py": k_py,
        "rel": rel,
        "max_abs_rel": np.nanmax(np.abs(rel)),
        "median_abs_rel": np.nanmedian(np.abs(rel)),
    }


# ============================================================
# DENSE TABLE GENERATION
# ============================================================

def generate_dense_table(
    logrho_min,
    logrho_max,
    nrho,
    T_min,
    T_max,
    ntemp,
    dust_coeff,
    gas_rho_grid,
    gas_T_grid,
    eG,
):
    logrho = np.linspace(
        logrho_min,
        logrho_max,
        nrho,
        dtype=np.float64,
    )

    logT = np.linspace(
        np.log10(T_min),
        np.log10(T_max),
        ntemp,
        dtype=np.float64,
    )

    rho_cgs = 10.0**logrho
    T_K = 10.0**logT

    kappa = np.empty(
        (nrho, ntemp),
        dtype=np.float64,
    )

    for ir, rho in enumerate(rho_cgs):
        for it, T in enumerate(T_K):
            kappa[ir, it] = semenov_rosseland(
                rho,
                T,
                dust_coeff,
                gas_rho_grid,
                gas_T_grid,
                eG,
            )

    if not np.all(np.isfinite(kappa)):
        bad = np.argwhere(~np.isfinite(kappa))[0]
        raise RuntimeError(
            "Non-finite opacity generated at "
            f"rho={rho_cgs[bad[0]]:.8e}, T={T_K[bad[1]]:.8e}."
        )

    if np.any(kappa <= 0.0):
        bad = np.argwhere(kappa <= 0.0)[0]
        raise RuntimeError(
            "Zero/negative opacity generated inside requested table domain at "
            f"rho={rho_cgs[bad[0]]:.8e}, T={T_K[bad[1]]:.8e}. "
            "Reduce the requested domain or inspect Semenov validity."
        )

    logkappa = np.log10(kappa)

    return (
        logrho,
        logT,
        rho_cgs,
        T_K,
        kappa,
        logkappa,
    )


def write_ascii_table(
    path,
    logrho,
    logT,
    logkappa,
):
    """
    Write a compact matrix-style ASCII table.

    Format:
        # SEMENOV...
        NRHO NTEMP
        <logrho values, one line>
        <logT values, one line>
        then NRHO rows, each containing NTEMP log10(kappa_R) values.

    This is intentionally easy to parse from C later.
    """
    path = Path(path)

    with path.open("w", encoding="utf-8") as f:
        f.write("# Semenov et al. (2003) Rosseland mean opacity\n")
        f.write("# dust model: nrm / h / s\n")
        f.write("# axes: log10(rho[g cm^-3]), log10(T[K])\n")
        f.write("# values: log10(kappa_R[cm^2 g^-1])\n")
        f.write(f"{logrho.size:d} {logT.size:d}\n")

        f.write(" ".join(f"{x:.16e}" for x in logrho))
        f.write("\n")

        f.write(" ".join(f"{x:.16e}" for x in logT))
        f.write("\n")

        for row in logkappa:
            f.write(" ".join(f"{x:.16e}" for x in row))
            f.write("\n")


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate a dense Semenov et al. (2003) nrm/h/s "
            "Rosseland opacity table using Python only."
        )
    )

    parser.add_argument(
        "--opacity-f",
        type=Path,
        default=Path("opacity.f"),
        help="Path to original Semenov opacity.f.",
    )

    parser.add_argument(
        "--gas-table",
        type=Path,
        default=Path("kR_h2001.dat"),
        help="Path to original kR_h2001.dat.",
    )

    parser.add_argument(
        "--reference",
        type=Path,
        default=None,
        help=(
            "Optional original FORTRAN kR.out for regression validation. "
            "The standard supplied kR.out corresponds to rho=1e-10 g cm^-3."
        ),
    )

    parser.add_argument(
        "--reference-rho",
        type=float,
        default=1.0e-10,
        help="Density used by the optional kR.out regression.",
    )

    parser.add_argument(
        "--logrho-min",
        type=float,
        default=-18.0,
        help="Minimum log10(rho[g cm^-3]) of dense output table.",
    )

    parser.add_argument(
        "--logrho-max",
        type=float,
        default=-7.0,
        help="Maximum log10(rho[g cm^-3]) of dense output table.",
    )

    parser.add_argument(
        "--nrho",
        type=int,
        default=241,
        help="Number of log-density grid points.",
    )

    parser.add_argument(
        "--T-min",
        type=float,
        default=5.0,
        help="Minimum temperature [K].",
    )

    parser.add_argument(
        "--T-max",
        type=float,
        default=10000.0,
        help="Maximum temperature [K].",
    )

    parser.add_argument(
        "--ntemp",
        type=int,
        default=801,
        help="Number of log-temperature grid points.",
    )

    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("semenov_rosseland_nrm_h_s"),
        help="Output path prefix.",
    )

    args = parser.parse_args()

    if args.logrho_min < np.log10(GAS_RHO_MIN):
        print(
            "WARNING: requested density range extends below the formal "
            "gop() limit. This may still be harmless where dust dominates, "
            "but a hot-gas query there can return zero."
        )

    if args.logrho_max > np.log10(GAS_RHO_MAX):
        raise ValueError(
            "Requested density maximum exceeds the original gop() maximum "
            f"{GAS_RHO_MAX:.3e} g cm^-3."
        )

    if args.T_min < 1.0:
        raise ValueError("Semenov COP() returns zero below 1 K.")

    if args.T_max > GAS_T_MAX:
        raise ValueError(
            "Requested T_max exceeds the original Semenov gas-table maximum "
            f"{GAS_T_MAX:.0f} K."
        )

    print("=" * 78)
    print("SEMENOV ET AL. (2003) ROSSeland TABLE GENERATOR")
    print("=" * 78)

    print()
    print("Reading original Semenov data:")
    print(" opacity.f      =", args.opacity_f)
    print(" kR_h2001.dat   =", args.gas_table)

    dust_coeff = parse_nrm_h_s_rosseland_coefficients(
        args.opacity_f
    )

    gas_rho_grid, gas_T_grid = parse_gop_grids(
        args.opacity_f
    )

    dummy, eG = read_kR_h2001(
        args.gas_table
    )

    print()
    print("Parsed model:")
    print(" dust model      = nrm / h / s")
    print(" opacity kind    = Rosseland")
    print(" dust fits       =", dust_coeff.shape)
    print(" gas grid        =", eG.shape)
    print(" gas dummy value =", dummy)
    print(
        " gas rho range  = "
        f"{gas_rho_grid[-1]:.6e} ... {gas_rho_grid[0]:.6e} g cm^-3"
    )
    print(
        " gas T range    = "
        f"{gas_T_grid[0]:.2f} ... {gas_T_grid[-1]:.2f} K"
    )

    if args.reference is not None:
        print()
        print("Regression against original FORTRAN kR.out:")

        check = validate_against_kR_out(
            args.reference,
            args.reference_rho,
            dust_coeff,
            gas_rho_grid,
            gas_T_grid,
            eG,
        )

        print(
            " max |relative difference|    = "
            f"{check['max_abs_rel']:.12e}"
        )
        print(
            " median |relative difference| = "
            f"{check['median_abs_rel']:.12e}"
        )
        print(
            " NOTE: kR.out is text-rounded, so ~1e-5 differences "
            "are expected even for a correct reproduction."
        )

    print()
    print("Generating dense table:")
    print(
        " log10 rho      = "
        f"[{args.logrho_min:.3f}, {args.logrho_max:.3f}]"
    )
    print(
        " T              = "
        f"[{args.T_min:.3f}, {args.T_max:.3f}] K"
    )
    print(" NRHO           =", args.nrho)
    print(" NTEMP          =", args.ntemp)

    (
        logrho,
        logT,
        rho_cgs,
        T_K,
        kappa,
        logkappa,
    ) = generate_dense_table(
        args.logrho_min,
        args.logrho_max,
        args.nrho,
        args.T_min,
        args.T_max,
        args.ntemp,
        dust_coeff,
        gas_rho_grid,
        gas_T_grid,
        eG,
    )

    print()
    print("Generated opacity range:")
    print(" kappa_R min =", np.min(kappa), "cm^2 g^-1")
    print(" kappa_R max =", np.max(kappa), "cm^2 g^-1")
    print(" finite      =", np.all(np.isfinite(kappa)))
    print(" positive    =", np.all(kappa > 0.0))

    npz_path = args.output_prefix.with_suffix(".npz")
    dat_path = args.output_prefix.with_suffix(".dat")
    fig_path = args.output_prefix.parent / (
        "fig_" + args.output_prefix.name + ".png"
    )

    np.savez_compressed(
        npz_path,
        logrho=logrho,
        logT=logT,
        rho_cgs=rho_cgs,
        T_K=T_K,
        kappa_cgs=kappa,
        logkappa=logkappa,
        model="nrm_h_s",
        opacity_kind="Rosseland",
    )

    write_ascii_table(
        dat_path,
        logrho,
        logT,
        logkappa,
    )

    # --------------------------------------------------------
    # Diagnostic map
    # --------------------------------------------------------

    fig, ax = plt.subplots(figsize=(8.2, 5.8))

    im = ax.pcolormesh(
        logT,
        logrho,
        logkappa,
        shading="auto",
    )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(
        r"$\log_{10}\kappa_R\ [{\rm cm^2\,g^{-1}}]$"
    )

    ax.set_xlabel(
        r"$\log_{10}T\ [{\rm K}]$"
    )

    ax.set_ylabel(
        r"$\log_{10}\rho\ [{\rm g\,cm^{-3}}]$"
    )

    ax.set_title(
        "Semenov et al. (2003): nrm / h / s Rosseland opacity"
    )

    plt.tight_layout()
    plt.savefig(
        fig_path,
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)

    print()
    print("Written:")
    print(" ", npz_path)
    print(" ", dat_path)
    print(" ", fig_path)

    print()
    print("=" * 78)
    print("TABLE GENERATION: PASS")
    print("=" * 78)


if __name__ == "__main__":
    main()
