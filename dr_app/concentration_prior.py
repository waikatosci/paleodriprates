"""
Lognormal concentration prior for stochastic aqueous chemistry
==============================================================

Companion to params.py — kept separate to avoid modifying the
original file. Import as:
    from concentration_prior import ConcentrationPrior

Both Ca_aq and each TE_aq use the same type. The prior is specified
in lognormal space (mu_ln, sigma_ln) and can be constructed from:
    - a monitoring time series  (from_series)
    - manual mean + SD          (from_mean_sd)
    - direct lognormal params   (constructor)

Sampling always returns ppb regardless of the input unit.
"""

from dataclasses import dataclass
import numpy as np


@dataclass
class ConcentrationPrior:
    """Lognormal prior for an aqueous concentration.

    Parameters
    ----------
    mu_ln : float
        Mean of ln(concentration).
    sigma_ln : float
        Std dev of ln(concentration).
    unit : str
        Unit that mu_ln was derived in: 'ppb' or 'ppm'.
        Sampling converts to ppb internally.
    source : str
        'csv' or 'manual' — for logging only.
    """
    mu_ln:    float
    sigma_ln: float
    unit:     str  = 'ppb'
    source:   str  = 'manual'

    def sample(self, n: int, rng=None) -> np.ndarray:
        """Return n samples in ppb."""
        rng = rng or np.random.default_rng()
        samples = rng.lognormal(self.mu_ln, self.sigma_ln, n)
        if self.unit == 'ppm':
            samples = samples * 1000.0
        return samples

    @property
    def mean_ppb(self) -> float:
        """Lognormal mean, converted to ppb."""
        mean = np.exp(self.mu_ln + 0.5 * self.sigma_ln ** 2)
        return mean * 1000.0 if self.unit == 'ppm' else mean

    @property
    def median_ppb(self) -> float:
        """Lognormal median, converted to ppb."""
        med = np.exp(self.mu_ln)
        return med * 1000.0 if self.unit == 'ppm' else med

    @property
    def cv(self) -> float:
        """Coefficient of variation (dimensionless)."""
        return float(np.sqrt(np.exp(self.sigma_ln ** 2) - 1))

    @classmethod
    def from_series(cls, values_ppb: np.ndarray,
                    source: str = 'csv') -> 'ConcentrationPrior':
        """Fit lognormal from a 1-D array of concentrations in ppb.

        Drops NaN, non-finite, and non-positive values before fitting.
        Raises ValueError if fewer than 3 valid values remain.
        """
        v = np.asarray(values_ppb, dtype=float)
        v = v[np.isfinite(v) & (v > 0)]
        if len(v) < 3:
            raise ValueError(
                f'Need >= 3 finite positive values to fit prior, got {len(v)}')
        ln_v = np.log(v)
        return cls(
            mu_ln=float(ln_v.mean()),
            sigma_ln=float(ln_v.std(ddof=1)),   # sample std dev
            unit='ppb',
            source=source,
        )

    @classmethod
    def from_mean_sd(cls, mean_ppb: float, sd_ppb: float,
                     source: str = 'manual') -> 'ConcentrationPrior':
        """Convert mean/SD on linear scale to lognormal parameters.

        Applies a 1% CV floor to avoid degenerate zero-variance priors.
        """
        if mean_ppb <= 0:
            raise ValueError('mean_ppb must be > 0')
        sd_ppb = max(sd_ppb, mean_ppb * 0.01)   # floor at 1% CV
        sigma2 = np.log(1 + (sd_ppb / mean_ppb) ** 2)
        mu = np.log(mean_ppb) - 0.5 * sigma2
        return cls(
            mu_ln=float(mu),
            sigma_ln=float(np.sqrt(sigma2)),
            unit='ppb',
            source=source,
        )

    def __repr__(self) -> str:
        return (f'ConcentrationPrior(mu_ln={self.mu_ln:.4f}, '
                f'sigma_ln={self.sigma_ln:.4f}, '
                f'mean={self.mean_ppb:.2f} ppb, '
                f'CV={self.cv * 100:.1f}%, '
                f'source={self.source!r})')


# ── Default scalar fallbacks (ppb) ──────────────────────────────────────
CA_AQ_DEFAULT_PPB = 62_000.0    # 62 ppm → ppb
NI_AQ_DEFAULT_PPB =      4.37
CO_AQ_DEFAULT_PPB =      0.30
CU_AQ_DEFAULT_PPB =      1.00

# ── Default number of concentration samples for stochastic mode ─────────
N_CONC_SAMPLES = 50
