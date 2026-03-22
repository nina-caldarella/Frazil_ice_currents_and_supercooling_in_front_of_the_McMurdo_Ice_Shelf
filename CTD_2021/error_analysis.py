import xarray as xr
import numpy as np
# ----------------------------
# EOS-80 freezing point (ITS-90)
# t_f(°C, ITS-90) = (a0*S + a1*S*sqrt(S) + a2*S^2 + b*p) × 0.99976
# a0..b are UNESCO / Millero coefficients (originally on IPTS-68).
# Multiply by 0.99976 to convert to ITS-90.  (1/1.00024 ≈ 0.99976)
# Refs: UNESCO 1983; Millero & Leung (1976); CSIRO sw_fp.m implementation.
# ----------------------------
from seawater import fp

# ----------------------------
# Practical Salinity from conductivity (PSS-78)
# Try python-seawater (EOS-80) first; fall back to TEOS-10 gsw.SP_from_C.
# Both implement the PSS-78 practical-salinity algorithm from UNESCO (1983).
# ----------------------------
from gsw import SP_from_C

## density
from seawater import pden

# ----------------------------
# Monte-Carlo uncertainty propagation for SP and EOS-80 freezing point
# ----------------------------
def propagate_uncertainty_eos80(
    C_Sm, T90_C, p_dbar,
    uC_Sm, uT_C, uP_dbar,
    n_samples=20000, seed=0, clip_C_positive=True, return_fd_check=True
):
    """
    Inputs may be scalars or arrays (they will be broadcast to common shape).
    Instrument 1-sigma uncertainties: uC (mS/cm), uT (°C ITS-90), uP (dbar).
    Returns dict with central values and 1-sigma uncertainties from Monte Carlo.
    """
    rng = np.random.default_rng(seed)

    # Broadcast to common shape
    C0 = np.asarray(C_Sm, float)
    T0 = np.asarray(T90_C, float)
    P0 = np.asarray(p_dbar, float)
    shape = np.broadcast_shapes(C0.shape, T0.shape, P0.shape)
    C0 = np.broadcast_to(C0, shape)
    T0 = np.broadcast_to(T0, shape)
    P0 = np.broadcast_to(P0, shape)

    # Draw samples
    C_s = rng.normal(C0, uC_Sm, size=(n_samples,) + shape)
    T_s = rng.normal(T0, uT_C,     size=(n_samples,) + shape)
    P_s = rng.normal(P0, uP_dbar,  size=(n_samples,) + shape)

    if clip_C_positive:
        C_s = np.maximum(C_s, np.finfo(float).tiny)  # avoid non-positive C on log/ratio ops

    # Convert C to units mS/cm for GSW toolbox
    C_s_mScm = 10 * C_s

    # Evaluate SP and freezing point for each draw
    # (loop by chunk to keep memory in check if your arrays are large)
    SP_s = np.empty_like(C_s)
    rho_s = np.empty_like(C_s)
    for k in range(n_samples):
        SP_s[k] = SP_from_C(C_s_mScm[k], T_s[k], P_s[k]) # replaced with original GSW package
        rho_s[k] = pden(SP_s[k], T_s[k]/ 0.99976, P_s[k])
        
    Tf_s = fp(SP_s, P_s) # replaced with original EOS80 package

    # Central values at nominal inputs
    SP0 = SP_from_C(10*C0, T0, P0)
    Tf0 = fp(SP0, P0)
    rho0 = pden(SP0, T0/ 0.99976, P0)

    # Monte-Carlo 1-sigma
    u_SP = SP_s.std(axis=0, ddof=1)
    u_rho = rho_s.std(axis=0, ddof=1)
    u_Tf = Tf_s.std(axis=0, ddof=1)
    
    depth_dummy = np.arange(len(C_Sm))

    ds = xr.Dataset({
        "SP": (("depth_dummy"), SP0),        # Practical Salinity (PSS-78)
        "rho": (("depth_dummy"), rho0),      # Potential density
        "Tf": (("depth_dummy"), Tf0),        # Freezing point (°C, ITS-90)
        "u_SP": (("depth_dummy"), u_SP),     # 1-sigma uncertainty in SP (PSS-78 units)
        "u_Tf": (("depth_dummy"), u_Tf),     # 1-sigma uncertainty in Tf (°C)
        "u_rho": (("depth_dummy"), u_rho),      # Potential density
         },
         coords = {"depth_dummy": depth_dummy,
         "n_samples": n_samples})

    return ds
