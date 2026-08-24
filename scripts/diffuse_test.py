import numpy as np
import fargopy as fp
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm



OUTPUT_DIR = "/home/alejo/research/fargo3d-FLD/outputs/rad_probe_00/"
SNAP = 2
# ============================================================
# SIMULATION AND UNITS
# ============================================================



sim = fp.Simulation(output_dir=OUTPUT_DIR)

sim.units("CGS")
sim.set_units(
    UM=1 * sim.MSUN,
    UL=5.2 * sim.AU
)


# ============================================================
# LOAD DATA
# ============================================================

profiles_E = []
profiles_u = []
times = []

r_grid = None

for snap in SNAPS:

    fields = sim.load_field(
        fields=["energyrad"],
        snapshot=snap,
        coords="spherical"
    )

    # FARGOpy meshes: always [0]
    phi_mesh   = np.asarray(fields.var1_mesh[0])
    r_mesh     = np.asarray(fields.var2_mesh[0])
    theta_mesh = np.asarray(fields.var3_mesh[0])

    erad = np.asarray(fields.energyrad_mesh[0])

    if erad.ndim != 3:
        raise RuntimeError(
            f"energyrad must be 3D, got shape={erad.shape}"
        )

    nz, ny, nx = erad.shape

    # ========================================================
    # RADIAL GRID
    # ========================================================

    if r_grid is None:

        print("\nMESH INFORMATION")
        print("-"*72)
        print(f"energyrad shape = {erad.shape}")
        print(f"var1_mesh shape = {phi_mesh.shape}")
        print(f"var2_mesh shape = {r_mesh.shape}")
        print(f"var3_mesh shape = {theta_mesh.shape}")

        # Scalar field ordering:
        # axis 0 -> theta
        # axis 1 -> r
        # axis 2 -> phi
        r_grid = np.asarray(
            r_mesh[0, :, 0],
            dtype=float
        )

        print("\nRADIAL GRID")
        print("-"*72)
        print(f"Nr       = {r_grid.size}")
        print(f"Expected = {ny}")
        print(f"r_min    = {r_grid.min():.12f}")
        print(f"r_max    = {r_grid.max():.12f}")
        print(f"r first  = {r_grid[:5]}")
        print(f"r last   = {r_grid[-5:]}")
        print(f"dr_min   = {np.min(np.diff(r_grid)):.12e}")
        print(f"dr_max   = {np.max(np.diff(r_grid)):.12e}")

        if r_grid.size != ny:
            raise RuntimeError(
                f"Radial grid has {r_grid.size} cells, "
                f"expected Ny={ny}."
            )

        if not np.all(np.diff(r_grid) > 0.0):
            raise RuntimeError(
                "Spherical radial coordinate is not monotonically increasing."
            )

    # ========================================================
    # ANGULAR AVERAGE
    # ========================================================

    Er_r = np.mean(
        erad,
        axis=(0, 2)
    )

    if Er_r.size != r_grid.size:
        raise RuntimeError(
            f"Radial EnergyRad profile has {Er_r.size} cells, "
            f"but radial grid has {r_grid.size}."
        )

    # ========================================================
    # u = r E_rad
    # ========================================================

    u = r_grid*Er_r

    profiles_E.append(Er_r.copy())
    profiles_u.append(u.copy())

    t = snap*DT_OUTPUT
    times.append(t)

    print(
        f"snapshot={snap:3d}  "
        f"t={t:10.5f}  "
        f"Er_min={np.min(erad):.8e}  "
        f"Er_max={np.max(erad):.8e}  "
        f"u_max={np.max(u):.8e}"
    )

profiles_E = np.asarray(profiles_E)
profiles_u = np.asarray(profiles_u)
times = np.asarray(times)
r_grid = np.asarray(r_grid)



# ============================================================
# RADIAL CELL WIDTHS
# ============================================================

r_edges = np.empty(r_grid.size + 1)

r_edges[1:-1] = 0.5 * (r_grid[:-1] + r_grid[1:])

r_edges[0] = (
    r_grid[0]
    - 0.5 * (r_grid[1] - r_grid[0])
)

r_edges[-1] = (
    r_grid[-1]
    + 0.5 * (r_grid[-1] - r_grid[-2])
)

dr = np.diff(r_edges)

# ============================================================
# GAUSSIAN MOMENTS OF u = r E_rad
# ============================================================

centers = []
variances = []
sigmas = []
u_integrals = []

for u in profiles_u:

    u_pos = np.maximum(u, 0.0)

    norm = np.sum(u_pos * dr)

    if norm <= 0.0:
        raise RuntimeError(
            "Integral of u = r E_rad is non-positive."
        )

    center = np.sum(
        r_grid * u_pos * dr
    ) / norm

    variance = np.sum(
        (r_grid - center)**2 * u_pos * dr
    ) / norm

    centers.append(center)
    variances.append(variance)
    sigmas.append(np.sqrt(variance))
    u_integrals.append(norm)

centers = np.asarray(centers)
variances = np.asarray(variances)
sigmas = np.asarray(sigmas)
u_integrals = np.asarray(u_integrals)

# ============================================================
# TOTAL RADIATION ENERGY
# ============================================================

# For spherical geometry:
#
#   E_tot ∝ ∫ E_rad r^2 dr
#
# Angular factors cancel in the relative error.

radiation_integrals = np.sum(
    profiles_E *
    r_grid[None, :]**2 *
    dr[None, :],
    axis=1
)

radiation_error = (
    radiation_integrals - radiation_integrals[0]
) / radiation_integrals[0]

# ============================================================
# ANALYTIC DIFFUSION LAW
# ============================================================

sigma_num0 = sigmas[0]
center_num0 = centers[0]

sigma_exact = np.sqrt(
    sigma_num0**2 + 2.0 * D * times
)

variance_growth_num = (
    variances - variances[0]
)

variance_growth_exact = (
    2.0 * D * times
)

# ============================================================
# MEASURE D
# ============================================================

slope, intercept = np.polyfit(
    times,
    variance_growth_num,
    1
)

D_measured = 0.5 * slope

relative_D_error = (
    D_measured / D - 1.0
)

# Instantaneous estimate from every snapshot

D_inst = np.full(
    times.shape,
    np.nan
)

D_inst[1:] = (
    variances[1:] - variances[0]
) / (
    2.0 * times[1:]
)

# ============================================================
# CONSERVATION DIAGNOSTICS
# ============================================================

center_shift = (
    centers - centers[0]
)

u_integral_error = (
    u_integrals - u_integrals[0]
) / u_integrals[0]

u_peak_num = np.max(
    profiles_u,
    axis=1
)

u_peak0 = u_peak_num[0]

u_peak_exact = (
    u_peak0 *
    sigma_num0 /
    sigma_exact
)

# ============================================================
# TERMINAL REPORT
# ============================================================

print("\n" + "="*76)
print("GAUSSIAN SPHERICAL DIFFUSION TEST")
print("="*76)

print(f"Input D                         = {D:.10e}")
print(f"Measured D                      = {D_measured:.10e}")
print(f"Relative D error                = {relative_D_error:+.6e}")

print()

print(f"Configured sigma0               = {SIGMA0:.10e}")
print(f"Measured sigma0                 = {sigma_num0:.10e}")

print()

print(f"Configured center               = {RC0:.10e}")
print(f"Measured center                 = {center_num0:.10e}")

print()

print(
    f"Maximum |center-center0|        = "
    f"{np.max(np.abs(center_shift)):.6e}"
)

print(
    f"Maximum radiation energy error = "
    f"{np.max(np.abs(radiation_error)):.6e}"
)

print(
    f"Maximum u-integral change       = "
    f"{np.max(np.abs(u_integral_error)):.6e}"
)

print("="*76)

print(
    f"\n{'snap':>5} "
    f"{'t':>10} "
    f"{'sigma_num':>14} "
    f"{'sigma_exact':>14} "
    f"{'D_inst':>14} "
    f"{'center':>12} "
    f"{'dE/E0':>14}"
)

for n, snap in enumerate(SNAPS):

    print(
        f"{snap:5d} "
        f"{times[n]:10.4f} "
        f"{sigmas[n]:14.7e} "
        f"{sigma_exact[n]:14.7e} "
        f"{D_inst[n]:14.7e} "
        f"{centers[n]:12.7f} "
        f"{radiation_error[n]:14.7e}"
    )

# ============================================================
# FIGURE
# ============================================================

fig, ax = plt.subplots(
    2, 2,
    figsize=(12.5, 9.0),
    constrained_layout=True
)

# ============================================================
# 1. NUMERICAL PROFILES
# ============================================================

for snap in SNAPS_TO_PLOT:

    if snap not in SNAPS:
        continue

    n = SNAPS.index(snap)

    ax[0,0].plot(
        r_grid,
        profiles_u[n],
        label=f"n={snap}"
    )

ax[0,0].set_xlabel(r"$r$")
ax[0,0].set_ylabel(r"$u=rE_{\rm rad}$")
ax[0,0].set_title("Numerical diffusion")
ax[0,0].legend()
ax[0,0].grid(alpha=0.25)

# ============================================================
# 2. NUMERICAL VS ANALYTIC
# ============================================================

for snap in SNAPS_TO_PLOT:

    if snap not in SNAPS:
        continue

    n = SNAPS.index(snap)

    sigma_t = sigma_exact[n]

    u_exact = (
        u_peak0 *
        sigma_num0 / sigma_t *
        np.exp(
            -0.5 *
            (r_grid - center_num0)**2 /
            sigma_t**2
        )
    )

    ax[0,1].plot(
        r_grid,
        profiles_u[n],
        linewidth=1.4,
        label=f"Numerical n={snap}"
    )

    ax[0,1].plot(
        r_grid,
        u_exact,
        "--",
        linewidth=1.0
    )

ax[0,1].set_xlabel(r"$r$")
ax[0,1].set_ylabel(r"$u=rE_{\rm rad}$")
ax[0,1].set_title("Numerical vs analytic")
ax[0,1].legend(fontsize=8)
ax[0,1].grid(alpha=0.25)

# ============================================================
# 3. VARIANCE GROWTH
# ============================================================

ax[1,0].plot(
    times,
    variance_growth_num,
    "o-",
    label="FARGO3D"
)

ax[1,0].plot(
    times,
    variance_growth_exact,
    "--",
    label=r"$2Dt$"
)

ax[1,0].set_xlabel(r"$t$")
ax[1,0].set_ylabel(
    r"$\sigma^2(t)-\sigma^2(0)$"
)

ax[1,0].set_title(
    rf"$D_{{meas}}={D_measured:.3e}$"
)

ax[1,0].legend()
ax[1,0].grid(alpha=0.25)

# ============================================================
# 4. CONSERVATION
# ============================================================

ax[1,1].plot(
    times,
    radiation_error,
    "o-",
    label=r"$\Delta E_{\rm rad}/E_{\rm rad,0}$"
)

ax[1,1].plot(
    times,
    center_shift,
    "s-",
    label=r"$r_c-r_{c,0}$"
)

ax[1,1].axhline(
    0.0,
    linestyle="--",
    linewidth=0.8
)

ax[1,1].set_xlabel(r"$t$")
ax[1,1].set_ylabel("Relative error / displacement")
ax[1,1].set_title("Conservation")
ax[1,1].legend()
ax[1,1].grid(alpha=0.25)

plt.show()