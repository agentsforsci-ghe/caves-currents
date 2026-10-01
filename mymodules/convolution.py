"""Convolve a forcing with impulse-response kernels along an AMOC pathway.

The response to a forcing F (decadal source-region d18O anomaly) with a
decadal kernel h is

    A(t) = sum_k c * F(t - k) * h(k)

c = SCALE = 1. The 2025.05 convolution notebook used c = 10, reasoning that
the dye pulse was 10 years at 1 normalised unit per year and F is a per-year
source anomaly. RESOLVED 2026-10-01 in favour of c = 1: a forcing held
constant for 500 yr reproduces the meltmodel paper's site anomalies with
c = 1 (median ratio 0.98 over 4 scenarios x 2 modes x 9 sites) and gives ten
times them with c = 10 (scripts/validate_constant_forcing.py). Laura's
decision. Earlier evidence pointed the same way. The 50 decadal pulse responses already add up to the constant-input
equilibrium: 0.85-1.16 per site for cold and zonal, and the notebook's own
constant runs agree. The source-region dye also decays smoothly from the
parent's pulse decade into the pulse run, with no 10x jump. With SCALE = 10,
a steady forcing therefore gives 10x the equilibrium site weights of Endres
et al. (2026b). SCALE = 1 reproduces those weights. Shares, dominance,
variance shares and the spider shapes do not depend on SCALE; anomaly values
and SDs scale with it. Pass `--scale` to scripts/run_convolution.py to change the pipeline itself.
A physical check: with SCALE = 10 the surface anomaly near the sources reaches
about -44 per mil, more negative than the -35 per mil meltwater end-member
itself, which no mixing can produce. With SCALE = 1 the maximum is about -4.4.

For a pathway that switches modes, the forcing of each decade goes into the
kernel of the mode active in that decade, and the responses are added. The
tail of a pulse released under one mode therefore keeps evolving after the
switch, which is what the notebook did by summing the segment convolutions.
"""

import numpy as np
from scipy.signal import fftconvolve

from myconfig.PATHWAYS import PATHWAYS

SCALE = 1.0   # no extra factor; validated against the meltmodel paper, see above


def pathway_modes(t_model, pathway):
    """Mode name for every model time (years, negative = BP)."""
    schedule = PATHWAYS[pathway] if isinstance(pathway, str) else pathway
    bp = -np.asarray(t_model)
    modes = np.empty(bp.shape, dtype=object)
    if len(schedule) == 1 and schedule[0][0] is None:
        modes[:] = schedule[0][2]
        return modes
    for older, younger, mode in schedule:
        modes[(bp <= older) & (bp > younger)] = mode
    modes[bp > schedule[0][0]] = schedule[0][2]        # spin-up before first segment
    modes[bp <= schedule[-1][1]] = schedule[-1][2]     # continue after last segment
    return modes


def convolve_series(forcing, point_kernel, pathway, scale=SCALE):
    """Regional contributions at one site.

    forcing       DataFrame (index model time, one column per region code)
    point_kernel  DataArray (mode, dye, lag), dyes in region-code order
    Returns an array (time, region) with the same time axis as `forcing`.
    """
    F = forcing.to_numpy(float)
    n = F.shape[0]
    modes = pathway_modes(forcing.index.to_numpy(), pathway)
    out = np.zeros_like(F)
    for mode in np.unique(modes):
        mask = (modes == mode)[:, None]
        h = point_kernel.sel(mode=mode).transpose("dye", "lag").to_numpy()
        for r in range(F.shape[1]):
            out[:, r] += np.convolve(scale * F[:, r] * mask[:, 0], h[r], mode="full")[:n]
    return out


def convolve_field(forcing, kernel_field, pathway, keep=None, scale=SCALE):
    """Total surface anomaly field (time, latitude, longitude).

    kernel_field  DataArray (mode, dye, lag, latitude, longitude)
    keep          optional boolean mask over forcing times to return
    Uses FFT convolution along time for all grid cells at once.
    """
    F = forcing.to_numpy(float)
    n = F.shape[0]
    modes = pathway_modes(forcing.index.to_numpy(), pathway)
    ny, nx = kernel_field.sizes["latitude"], kernel_field.sizes["longitude"]
    total = np.zeros((n, ny, nx))
    for mode in np.unique(modes):
        mask = (modes == mode)
        k = np.nan_to_num(kernel_field.sel(mode=mode).transpose("dye", "lag", "latitude", "longitude").to_numpy())
        for r in range(F.shape[1]):
            x = (scale * F[:, r] * mask)[:, None, None]
            total += fftconvolve(x, k[r], mode="full", axes=0)[:n]
    land = np.isnan(kernel_field.isel(mode=0, dye=0, lag=0).to_numpy())
    total[:, land] = np.nan
    return total if keep is None else total[keep]
