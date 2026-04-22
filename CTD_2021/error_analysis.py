import xarray as xr
import numpy as np

# ----------------------------
# Monte-Carlo uncertainty propagation for SP and EOS-80 freezing point
# ----------------------------
from seawater import fp, pden

def pressure_uncertainty_eos80(SP, T90, p, uP, dp=0.01):
    """
    Depth-dependent absolute uncertainty from pressure only (EOS-80).
    """

    # Reference
    Tf0  = fp(SP, p)                          # ITS-90
    rho0 = pden(SP, T90 / 0.99976, p)         # pden expects IPTS-68

    # Numerical sensitivities
    dTf_dp = (fp(SP, p + dp) - fp(SP, p - dp)) / (2 * dp)
    drho_dp = (
        pden(SP, T90 / 0.99976, p + dp)
      - pden(SP, T90 / 0.99976, p - dp)
    ) / (2 * dp)

    # Absolute uncertainties
    uTf_p  = np.abs(dTf_dp)  * uP
    urho_p = np.abs(drho_dp) * uP

    return uTf_p, urho_p

from gsw import SP_from_C

def scalar_relative_uncertainty_eos80(
    C_Sm, T90, p,
    uC_Sm, uT90,
    n_samples=5000, seed=0, clip_C_positive=True
):
    """
    Scalar relative uncertainty in Tf and rho from conductivity and temperature (EOS-80).
    """

    rng = np.random.default_rng(seed)

    # Representative conditions
    C0  = float(np.nanmedian(C_Sm))
    T0  = float(np.nanmedian(T90))
    P0  = float(np.nanmedian(p))

    uC  = float(uC_Sm)
    uT  = float(uT90)

    # Monte Carlo
    C_s = rng.normal(C0, uC, size=n_samples)
    T_s = rng.normal(T0, uT, size=n_samples)

    if clip_C_positive:
        C_s = np.maximum(C_s, np.finfo(float).tiny)

    SP_s = SP_from_C(10.0 * C_s, T0, P0)

    Tf_s  = fp(SP_s, P0)
    rho_s = pden(SP_s, T_s / 0.99976, P0)

    SP0   = SP_from_C(10.0 * C0, T0, P0)
    Tf0   = fp(SP0, P0)
    rho0  = pden(SP0, T0 / 0.99976, P0)

    rTf_CT  = Tf_s.std(ddof=1)  / abs(Tf0)
    rrho_CT = rho_s.std(ddof=1) / abs(rho0)
    ruSP    = SP_s.std(ddof=1)  / abs(SP0)

    return rTf_CT, rrho_CT, ruSP 

def total_uncertainty_EOS80(
    SP, T90, p, uP,
    C_Sm, uC_Sm,
    uT90
):
    """
    Optimised EOS-80 uncertainty estimation matching TEOS-10 structure.
    """

    # Scalar relative uncertainty from C and T
    rTf_CT, rrho_CT, ruSP = scalar_relative_uncertainty_eos80(
        C_Sm, T90, p, uC_Sm, uT90
    )

    # Pressure contribution (absolute, depth-dependent)
    uTf_p, urho_p = pressure_uncertainty_eos80(
        SP, T90, p, uP
    )

    # Reference profiles
    Tf0  = fp(SP, p)
    rho0 = pden(SP, T90 / 0.99976, p)

    # Convert scalar relative → absolute
    uTf_CT  = rTf_CT  * np.abs(Tf0)
    urho_CT = rrho_CT * np.abs(rho0)

    # Quadrature sum (absolute, depth-dependent)
    uTf_total  = np.sqrt(uTf_p**2  + uTf_CT**2)
    urho_total = np.sqrt(urho_p**2 + urho_CT**2)

    return {
        "uTf":  uTf_total,   # absolute (depth)
        "urho": urho_total,  # absolute (depth)
        "ruSP": ruSP         # scalar relative
    }
    
# ----------------------------
# Monte-Carlo uncertainty propagation for TEOS-10
# ----------------------------
import gsw

def pressure_uncertainty_teos10(SA, CT, p, uP, dp=0.01):
    """
    Depth-dependent absolute uncertainty from pressure only.
    """

    # Reference
    Tf0  = gsw.CT_freezing(SA, p, 0.0)
    rho0 = gsw.rho(SA, CT, p)

    # Numerical sensitivities
    dTf_dp = (
        gsw.CT_freezing(SA, p + dp, 0.0)
        - gsw.CT_freezing(SA, p - dp, 0.0)
    ) / (2 * dp)

    drho_dp = (
        gsw.rho(SA, CT, p + dp)
        - gsw.rho(SA, CT, p - dp)
    ) / (2 * dp)

    # Absolute uncertainties
    uTf_p  = np.abs(dTf_dp)  * uP
    urho_p = np.abs(drho_dp) * uP

    return uTf_p, urho_p

def uncertainty_Tf_rho_SA(
    SA, CT, C_Sm, p, uCT, uC_Sm, lon, lat,
    n_samples=5000, seed=0, clip_C_positive=True
):
    """
    Scalar relative uncertainty from CT and SA.
    """

    rng = np.random.default_rng(seed)

    CT0 = float(np.nanmedian(CT))
    P0  = float(np.nanmedian(p))
    C0  = float(np.nanmedian(C_Sm))
    SA0 = float(np.nanmedian(SA))

    uCT = float(uCT)
        
    uC  = float(uC_Sm)

    # Monte Carlo on conductivity
    C_s = rng.normal(C0, uC, size=n_samples)
    if clip_C_positive:
        C_s = np.maximum(C_s, np.finfo(float).tiny)

    # C → SP → SA
    SP_s = gsw.SP_from_C(10.0 * C_s, CT0, P0)
    SA_s = gsw.SA_from_SP(SP_s, P0, lon, lat)

    Tf_s  = gsw.CT_freezing(SA_s, P0, 0.0)
    rho_s = gsw.rho(SA_s, CT0, P0)

    Tf0  = gsw.CT_freezing(SA0, P0, 0.0)
    rho0 = gsw.rho(SA0, CT0, P0)

    rTf_CTSA  = Tf_s.std(ddof=1)  / abs(Tf0)
    rrho_CTSA = rho_s.std(ddof=1) / abs(rho0)
    ruSA = SA_s.std(ddof=1) / abs(SA0)

    return rTf_CTSA, rrho_CTSA, ruSA
    
def total_uncertainty_TEOS10(
    SA, CT, p, uP,
    C_Sm, uC_Sm,
    uCT,
    lon, lat
):

    # --- Scalar relative CT+SA uncertainty ---
    rTf_CTSA, rrho_CTSA, ruSA = uncertainty_Tf_rho_SA(SA, CT, C_Sm, p, uCT, uC_Sm, lon, lat,
    )

    # --- Pressure contribution ---
    uTf_p, urho_p = pressure_uncertainty_teos10(
        SA, CT, p, uP
    )

    # Reference values
    Tf0  = gsw.CT_freezing(SA, p, 0.0)
    rho0 = gsw.rho(SA, CT, p)

    # Convert scalar relative → absolute
    uTf_CTSA  = rTf_CTSA  * np.abs(Tf0)
    urho_CTSA = rrho_CTSA * np.abs(rho0)

    # Final absolute uncertainty vs depth
    uTf_total  = np.sqrt(uTf_p**2  + uTf_CTSA**2)
    urho_total = np.sqrt(urho_p**2 + urho_CTSA**2)

    return {
        "uTf": uTf_total,            # absolute, depth-dependent
        "urho": urho_total,          # absolute, depth-dependent
        "ruSA": ruSA                 # scalar absolute
    }
