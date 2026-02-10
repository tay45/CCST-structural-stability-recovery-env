#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ccst_env_recovery_sweep_vnext.py

CCST toy visualization (single-file, "final" extended version)

Goals
-----
Produce paper-ready, reproducible toy plots that visualize:
  (1) baseline environment coupling (drift backbone)
  (2) recovery geometry (quad/asym)
  (3) environmental fluctuations (white / OU), acting on:
        - lambda (effective coupling rate), or
        - Phi (environmental drive term)

Diagnostics
-----------
Structural survival fraction:
  survival iff min_{t in [0,T]} eta(t) >= eta_crit
  R_struct(eta_crit) = fraction surviving across N trajectories

Outputs
-------
(A) Overlay plot (3-stage):
    (1) env-only (no recovery, no fluctuations)
    (2) + recovery (no fluctuations)
    (3) + recovery + fluctuations
(B) Grid plot:
    columns = env_sigma values
    rows    = recovery strength values (k for quad; scale for asym)
    with shared y-axis scale across all panels (paper-ready)

CLI
---
- You can run a single plot with --overlay or --grid (or both).
- Or use --preset <name> to run paper-meaningful configurations.
- Or use --run-all-presets to generate "one of each kind" (final test set).

Important modeling note (honesty)
---------------------------------
This is a toy model: its job is visual regime illustration, not micro-branching.
We keep the "drift vs recovery" separation explicit, and introduce fluctuations
as environment variability acting on a chosen target (lambda or Phi).

Model
-----
State: eta(t) in [0,1] (coherence/stability proxy)

Baseline deterministic backbone (per branch b):
  d eta/dt = -lambda_b * eta - gamma * Phi(eta)
  Phi(eta) = Phi0 + alpha*(1 - eta)

Recovery (optional):
  + k_side(eta) * (eta_star - eta)
  where:
    none: k_side=0
    quad: k_side = k
    asym: k_side = k_minus if eta<eta_star else k_plus

Environmental fluctuations (optional)
-------------------------------------
Targets:
  target=lambda:
    lambda_b -> lambda_b + delta_lambda(t)
    enters drift as:  -(delta_lambda(t) * eta)

  target=Phi:
    Phi -> Phi + delta_Phi(t)
    enters drift as:  -(gamma * delta_Phi(t))

Fluctuation models:
  none:
    delta_* = 0
  white:
    EM increment:
      target=lambda: d eta += -(env_sigma * eta) dW
      target=Phi   : d eta += -(gamma * env_sigma) dW
    env_sigma is a diffusion amplitude (per sqrt(time))
  ou:
    OU state X(t) with correlation time tau and stationary std env_sigma:
      dX = -(1/tau)X dt + sqrt(2/tau)*env_sigma dW
    Insert X into drift as delta_lambda=X or delta_Phi=X.

Reproducibility
---------------
- Default uses CRN (common random numbers) across conditions to stabilize
  *comparisons*. Disable with --no-crn.

Examples
--------
# Overlay (default preset-like settings if you specify args)
python3 ccst_env_recovery_sweep_vnext.py --overlay --outdir figs

# Grid
python3 ccst_env_recovery_sweep_vnext.py --grid --env-sigma-list 0,0.1,0.2,0.3 --k-list 0,6,12,18 --outdir figs

# Presets
python3 ccst_env_recovery_sweep_vnext.py --preset grid_ou_phi_asym --outdir figs
python3 ccst_env_recovery_sweep_vnext.py --run-all-presets --outdir figs
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Literal, Tuple

import numpy as np
import matplotlib.pyplot as plt


RecoveryType = Literal["none", "quad", "asym"]
EnvModel = Literal["none", "white", "ou"]
FluctTarget = Literal["lambda", "Phi"]


# ============================================================
# Config
# ============================================================

@dataclass(frozen=True)
class Config:
    # Ensemble + time
    N: int = 20000
    T: float = 1.0
    dt: float = 1e-3
    eta0: float = 0.95

    # Baseline environment coupling distribution for lambda_b
    lambda_min: float = 0.02
    lambda_max: float = 0.55

    # Environment drive Phi(eta) = Phi0 + alpha*(1-eta)
    gamma: float = 0.8
    Phi0: float = 0.1
    alpha: float = 0.6

    # Recovery basin parameters
    eta_star: float = 0.85
    k: float = 12.0
    k_minus: float = 18.0
    k_plus: float = 6.0

    # Fluctuations
    env_sigma: float = 0.0
    tau: float = 0.05  # OU correlation time (if env_model=ou)

    # RNG controls
    seed: int = 12345
    crn: bool = True   # common random numbers across conditions by default


# ============================================================
# Helpers
# ============================================================

def _steps(T: float, dt: float) -> int:
    s = int(round(T / dt))
    if s <= 0:
        raise ValueError("T/dt must be positive.")
    return s


def _parse_csv_floats(s: str) -> List[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def _eta_crit_default(n: int = 10, lo: float = 0.50, hi: float = 0.95) -> List[float]:
    return [float(x) for x in np.linspace(lo, hi, n)]


def _seed_mix(base: int, tag: int) -> int:
    # deterministic mixing (avoid Python hash randomization)
    x = (base ^ (tag * 0x9E3779B1)) & 0xFFFFFFFF
    x ^= (x >> 16)
    x = (x * 0x85EBCA6B) & 0xFFFFFFFF
    x ^= (x >> 13)
    x = (x * 0xC2B2AE35) & 0xFFFFFFFF
    x ^= (x >> 16)
    return int(x)


def _global_ylim_from_panels(panels: List[np.ndarray]) -> Tuple[float, float]:
    # Panels are in percent scale.
    if not panels:
        return 0.0, 100.0
    a = np.vstack(panels)
    y_min = float(np.nanmin(a))
    y_max = float(np.nanmax(a))
    span = max(1e-6, y_max - y_min)
    y0 = max(0.0, y_min - 0.10 * span)
    y1 = min(100.0, y_max + 0.10 * span)
    if y1 - y0 < 1.0:
        mid = 0.5 * (y0 + y1)
        y0 = max(0.0, mid - 0.5)
        y1 = min(100.0, mid + 0.5)
    return y0, y1


# ============================================================
# Model pieces
# ============================================================

def phi_det(eta: np.ndarray, cfg: Config) -> np.ndarray:
    return cfg.Phi0 + cfg.alpha * (1.0 - eta)


def recovery_k_side(eta: np.ndarray, cfg: Config, recovery: RecoveryType) -> np.ndarray:
    if recovery == "none":
        return np.zeros_like(eta)
    if recovery == "quad":
        return np.full_like(eta, cfg.k, dtype=float)
    if recovery == "asym":
        out = np.empty_like(eta)
        m = eta < cfg.eta_star
        out[m] = cfg.k_minus
        out[~m] = cfg.k_plus
        return out
    raise ValueError(f"Unknown recovery={recovery}")


# ============================================================
# Simulation
# ============================================================

def simulate_eta_min(
    *,
    cfg: Config,
    recovery: RecoveryType,
    env_model: EnvModel,
    fluct_target: FluctTarget,
    env_sigma: float,
    tau: float,
    rng_seed: int,
) -> np.ndarray:
    """
    Vectorized simulation for N trajectories; returns eta_min on [0,T].

    White noise (Euler–Maruyama increment):
      target=lambda: eta <- eta - env_sigma*eta*sqrt(dt)*Z
      target=Phi   : eta <- eta - gamma*env_sigma*sqrt(dt)*Z

    OU noise:
      OU state X has stationary std env_sigma; correlation tau; exact discretization:
        X_{t+dt} = a*X_t + std_step*Z,
        a=exp(-dt/tau), std_step=env_sigma*sqrt(1-a^2)
      Insert:
        target=lambda: drift += -(X*eta)
        target=Phi   : drift += -(gamma*X)
    """
    if cfg.N <= 0:
        raise ValueError("N must be > 0")
    if cfg.dt <= 0 or cfg.T <= 0:
        raise ValueError("dt and T must be > 0")
    if env_sigma < 0:
        raise ValueError("env_sigma must be >= 0")
    if env_model == "ou" and tau <= 0:
        raise ValueError("tau must be > 0 for OU")

    steps = _steps(cfg.T, cfg.dt)
    rng = np.random.default_rng(rng_seed)

    lambdas = rng.uniform(cfg.lambda_min, cfg.lambda_max, size=cfg.N)
    eta = np.full(cfg.N, cfg.eta0, dtype=float)
    eta_min = eta.copy()

    # OU parameters
    if env_model == "ou" and env_sigma > 0:
        X = np.zeros(cfg.N, dtype=float)
        a = np.exp(-cfg.dt / tau)
        std_step = env_sigma * np.sqrt(max(0.0, 1.0 - a * a))
    else:
        X = None
        a = 0.0
        std_step = 0.0

    sqrt_dt = np.sqrt(cfg.dt)

    for _ in range(steps):
        # baseline drift
        Phi = phi_det(eta, cfg)
        drift = -(lambdas * eta) - (cfg.gamma * Phi)

        # recovery drift
        k_side = recovery_k_side(eta, cfg, recovery)
        drift = drift + k_side * (cfg.eta_star - eta)

        # OU update + insertion into drift
        if env_model == "ou" and env_sigma > 0:
            Z = rng.normal(0.0, 1.0, size=cfg.N)
            X = a * X + std_step * Z
            if fluct_target == "lambda":
                drift = drift - (X * eta)
            else:
                drift = drift - (cfg.gamma * X)

        # Euler step
        eta = eta + cfg.dt * drift

        # White noise increment (EM)
        if env_model == "white" and env_sigma > 0:
            Z = rng.normal(0.0, 1.0, size=cfg.N)
            if fluct_target == "lambda":
                eta = eta - (env_sigma * eta) * sqrt_dt * Z
            else:
                eta = eta - (cfg.gamma * env_sigma) * sqrt_dt * Z

        # bounds
        eta = np.clip(eta, 0.0, 1.0)
        eta_min = np.minimum(eta_min, eta)

    return eta_min


def rstruct_curve(
    *,
    cfg: Config,
    recovery: RecoveryType,
    env_model: EnvModel,
    fluct_target: FluctTarget,
    env_sigma: float,
    tau: float,
    eta_crit_list: List[float],
    seed_tag: int,
) -> List[float]:
    seed = cfg.seed if cfg.crn else _seed_mix(cfg.seed, seed_tag)
    eta_min = simulate_eta_min(
        cfg=cfg,
        recovery=recovery,
        env_model=env_model,
        fluct_target=fluct_target,
        env_sigma=env_sigma,
        tau=tau,
        rng_seed=seed,
    )
    return [float(np.mean(eta_min >= e)) for e in eta_crit_list]


# ============================================================
# Plotting
# ============================================================

def plot_overlay_three_stage(
    *,
    cfg: Config,
    eta_crit_list: List[float],
    recovery: RecoveryType,
    env_model: EnvModel,
    fluct_target: FluctTarget,
    env_sigma: float,
    tau: float,
    outpath: Path,
    paper_ready: bool = True,
    show_textbox: bool = True,
    textbox_mode: str = "compact",  # "compact" | "full"
) -> None:
    x = np.array(eta_crit_list, dtype=float)

    # Stage 1: env-only (no recovery, no fluctuations)
    y_env_only = np.array(
        rstruct_curve(
            cfg=cfg,
            recovery="none",
            env_model="none",
            fluct_target=fluct_target,
            env_sigma=0.0,
            tau=tau,
            eta_crit_list=eta_crit_list,
            seed_tag=1,
        ),
        dtype=float,
    )

    # Stage 2: + recovery (no fluctuations)
    y_recovery = np.array(
        rstruct_curve(
            cfg=cfg,
            recovery=recovery,
            env_model="none",
            fluct_target=fluct_target,
            env_sigma=0.0,
            tau=tau,
            eta_crit_list=eta_crit_list,
            seed_tag=2,
        ),
        dtype=float,
    )

    # Stage 3: + recovery + fluctuations
    y_recovery_fluct = np.array(
        rstruct_curve(
            cfg=cfg,
            recovery=recovery,
            env_model=env_model,
            fluct_target=fluct_target,
            env_sigma=env_sigma,
            tau=tau,
            eta_crit_list=eta_crit_list,
            seed_tag=3,
        ),
        dtype=float,
    )

    # styling
    if paper_ready:
        figsize = (7.2, 3.8)
        dpi = 240
        lw = 2.0
        ms = 3.8
        fs = 10
        fs_title = 11
    else:
        figsize = (10, 5)
        dpi = 160
        lw = 2.4
        ms = 4.8
        fs = 11
        fs_title = 13

    plt.figure(figsize=figsize, dpi=dpi)

    plt.plot(
        x,
        100 * y_env_only,
        marker="o",
        linewidth=lw,
        markersize=ms,
        label="(1) env only (no recovery, no fluctuations)",
    )
    plt.plot(
        x,
        100 * y_recovery,
        marker="o",
        linewidth=lw,
        markersize=ms,
        label=f"(2) + recovery ({recovery}), no fluctuations",
    )
    plt.plot(
        x,
        100 * y_recovery_fluct,
        marker="o",
        linewidth=lw,
        markersize=ms,
        label=f"(3) + fluctuations ({env_model}, target={fluct_target}, σ={env_sigma:g})",
    )

    y0, y1 = _global_ylim_from_panels(
        [100 * y_env_only, 100 * y_recovery, 100 * y_recovery_fluct]
    )
    plt.ylim(y0, y1)
    plt.xlim(float(x.min()) - 0.01, float(x.max()) + 0.01)
    plt.grid(True, alpha=0.25)
    plt.xlabel(r"Structural threshold $\eta_{\mathrm{crit}}$", fontsize=fs)
    plt.ylabel(r"$R_{\mathrm{struct}}$ (\%)", fontsize=fs)
    plt.title(
        "Three-stage overlay: baseline → recovery → recovery+fluctuations",
        fontsize=fs_title,
    )

    # parameter textbox (compact + unobtrusive)
    if show_textbox:
        if textbox_mode == "full":
            box = (
                rf"$N={cfg.N}$, $T={cfg.T}$, $\Delta t={cfg.dt}$, $\eta_0={cfg.eta0}$" + "\n"
                rf"$\lambda\sim U({cfg.lambda_min},{cfg.lambda_max})$, $\gamma={cfg.gamma}$, $\Phi_0={cfg.Phi0}$, $\alpha={cfg.alpha}$" + "\n"
                rf"$\eta^\ast={cfg.eta_star}$, quad $k={cfg.k}$, asym $(k_-,k_+)=( {cfg.k_minus},{cfg.k_plus})$" + "\n"
                rf"env={env_model}, target={fluct_target}, $\sigma={env_sigma:g}$, OU $\tau={tau}$; CRN={cfg.crn}, seed={cfg.seed}"
            )
        else:
            # compact (recommended for paper)
            box = (
                rf"$N={cfg.N}$, $T={cfg.T}$, $\Delta t={cfg.dt}$" + "\n"
                rf"$\lambda\sim U({cfg.lambda_min},{cfg.lambda_max})$, $\gamma={cfg.gamma}$" + "\n"
                rf"recovery={recovery}, env={env_model}, target={fluct_target}, $\sigma={env_sigma:g}$, $\tau={tau}$"
            )

        ax = plt.gca()
        ax.text(
            0.985,
            0.02,
            box,
            transform=ax.transAxes,
            fontsize=max(7, fs - 3),
            va="bottom",
            ha="right",
            bbox=dict(
                boxstyle="round,pad=0.18",
                facecolor="white",
                alpha=0.72,
                edgecolor="0.75",
            ),
        )

    plt.legend(frameon=True, fontsize=fs - 1, loc="best")
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath, bbox_inches="tight")
    plt.close()


def plot_grid_recovery_vs_env(
    *,
    cfg: Config,
    eta_crit_list: List[float],
    recovery: RecoveryType,
    env_model: EnvModel,
    fluct_target: FluctTarget,
    env_sigma_list: List[float],
    k_list: List[float],
    tau: float,
    outpath: Path,
    paper_ready: bool = True,
    annotate: Literal["none", "int", "2dp"] = "none",
    remove_panel_legends: bool = True,
    big_outer_labels: bool = True,
) -> None:
    """
    Paper-ready grid:
      columns = env_sigma
      rows    = recovery strength

    For quad: row uses k values.
    For asym: row uses a "scale factor" s; internally (k_minus,k_plus) := s*(base k_minus,k_plus).
    """
    x = np.array(eta_crit_list, dtype=float)
    env_sigma_list = [float(v) for v in env_sigma_list]
    k_list = [float(v) for v in k_list]

    nrows = len(k_list)
    ncols = len(env_sigma_list)

    # precompute all curves -> enforce shared y-limits
    panels: Dict[Tuple[int, int], np.ndarray] = {}
    all_panels: List[np.ndarray] = []

    base_km, base_kp = cfg.k_minus, cfg.k_plus

    for i, kval in enumerate(k_list):
        if recovery == "none":
            cfg_row = cfg
        elif recovery == "quad":
            cfg_row = Config(**{**asdict(cfg), "k": float(kval)})
        else:
            # asym: interpret kval as multiplier
            cfg_row = Config(
                **{
                    **asdict(cfg),
                    "k_minus": float(kval) * base_km,
                    "k_plus": float(kval) * base_kp,
                }
            )

        for j, es in enumerate(env_sigma_list):
            em = env_model if es > 0 else "none"
            y = np.array(
                rstruct_curve(
                    cfg=cfg_row,
                    recovery=recovery,
                    env_model=em,
                    fluct_target=fluct_target,
                    env_sigma=float(es),
                    tau=tau,
                    eta_crit_list=eta_crit_list,
                    seed_tag=1000 + 37 * i + j,
                ),
                dtype=float,
            )
            y_pct = 100.0 * y
            panels[(i, j)] = y_pct
            all_panels.append(y_pct)

    y0, y1 = _global_ylim_from_panels(all_panels)

    # styling
    if paper_ready:
        figsize = (2.45 * ncols + 1.25, 2.15 * nrows + 1.15)
        dpi = 320
        lw = 1.6
        ms = 3.0
        fs_tick = 8.5
        fs_title = 10.5
        fs_outer = 12
        fs_rowcol = 9.5
        fs_anno = 7.5
    else:
        figsize = (3.6 * ncols, 3.0 * nrows)
        dpi = 180
        lw = 2.0
        ms = 4.2
        fs_tick = 10
        fs_title = 12
        fs_outer = 14
        fs_rowcol = 11
        fs_anno = 9

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=figsize,
        dpi=dpi,
        sharex=True,
        sharey=True,
        squeeze=False,
        constrained_layout=True,
    )

    for i, kval in enumerate(k_list):
        for j, es in enumerate(env_sigma_list):
            ax = axes[i][j]
            y_pct = panels[(i, j)]

            ax.plot(x, y_pct, marker="o", linewidth=lw, markersize=ms)
            ax.grid(True, alpha=0.22)
            ax.set_xlim(float(x.min()) - 0.01, float(x.max()) + 0.01)
            ax.set_ylim(y0, y1)
            ax.tick_params(labelsize=fs_tick)

            # annotate if requested
            if annotate != "none":
                for xi, yi in zip(x, y_pct):
                    if annotate == "int":
                        lab = f"{int(round(yi))}"
                    else:
                        lab = f"{yi:.2f}"
                    ax.text(
                        xi,
                        yi,
                        lab,
                        ha="center",
                        va="bottom",
                        fontsize=fs_anno,
                        clip_on=True,
                    )

            # column headers (top row): env sigma
            if i == 0:
                ax.set_title(rf"$\sigma_{{\mathrm{{env}}}}={es:g}$", fontsize=fs_rowcol)

            # row labels (left side): recovery strength
            if j == 0:
                if recovery == "quad":
                    rowlab = rf"$k={kval:g}$"
                elif recovery == "asym":
                    rowlab = rf"$\times {kval:g}$"
                else:
                    rowlab = "none"
                ax.text(
                    -0.22,
                    0.50,
                    rowlab,
                    transform=ax.transAxes,
                    rotation=90,
                    va="center",
                    ha="center",
                    fontsize=fs_rowcol,
                )

            # outer axis labels only
            if i == nrows - 1:
                ax.set_xlabel(r"$\eta_{\mathrm{crit}}$", fontsize=fs_title)
            if j == 0:
                ax.set_ylabel(r"$R_{\mathrm{struct}}$ (\%)$", fontsize=fs_title)

            if not remove_panel_legends:
                ax.legend(fontsize=fs_tick)

    # big outer labels (paper-like)
    if big_outer_labels:
        fig.suptitle(
            rf"Grid: recovery={recovery} | env={env_model}, target={fluct_target} (OU $\tau={tau}$)",
            fontsize=fs_title + 1,
        )
        fig.text(
            0.5,
            0.01,
            "Env fluctuation strength  (columns: increasing $\\sigma_{\\mathrm{env}}$)",
            ha="center",
            va="bottom",
            fontsize=fs_outer,
        )
        fig.text(
            0.01,
            0.5,
            "Recovery strength  (rows: increasing basin strength)",
            ha="left",
            va="center",
            rotation=90,
            fontsize=fs_outer,
        )
    else:
        fig.suptitle(
            rf"Grid: recovery={recovery} | env={env_model}, target={fluct_target} (OU tau={tau})",
            fontsize=fs_title + 1,
        )

    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# Presets
# ============================================================

def preset_table() -> Dict[str, Dict]:
    """
    Each preset returns a dict specifying:
      mode: "overlay" or "grid"
      and corresponding argument overrides.
    """
    return {
        # --- Overlays (3-stage) ---
        "overlay_white_lambda_quad": dict(
            mode="overlay",
            recovery="quad",
            env_model="white",
            fluct_target="lambda",
            env_sigma=0.25,
            tau=0.05,
            paper_ready=True,
        ),
        "overlay_white_lambda_asym": dict(
            mode="overlay",
            recovery="asym",
            env_model="white",
            fluct_target="lambda",
            env_sigma=0.25,
            tau=0.05,
            paper_ready=True,
        ),
        "overlay_ou_phi_quad": dict(
            mode="overlay",
            recovery="quad",
            env_model="ou",
            fluct_target="Phi",
            env_sigma=0.25,
            tau=0.05,
            paper_ready=True,
        ),
        "overlay_ou_phi_asym": dict(
            mode="overlay",
            recovery="asym",
            env_model="ou",
            fluct_target="Phi",
            env_sigma=0.25,
            tau=0.05,
            paper_ready=True,
        ),

        # --- Grids (paper-ready) ---
        "grid_white_lambda_quad": dict(
            mode="grid",
            recovery="quad",
            env_model="white",
            fluct_target="lambda",
            env_sigma_list=[0.0, 0.10, 0.20, 0.30],
            k_list=[0.0, 6.0, 12.0, 18.0],
            tau=0.05,
            paper_ready=True,
        ),
        "grid_white_lambda_asym": dict(
            mode="grid",
            recovery="asym",
            env_model="white",
            fluct_target="lambda",
            env_sigma_list=[0.0, 0.10, 0.20, 0.30],
            k_list=[0.0, 0.5, 1.0, 1.5],  # multiplier on (k_minus,k_plus)
            tau=0.05,
            paper_ready=True,
        ),
        "grid_ou_phi_quad": dict(
            mode="grid",
            recovery="quad",
            env_model="ou",
            fluct_target="Phi",
            env_sigma_list=[0.0, 0.10, 0.20, 0.30],
            k_list=[0.0, 6.0, 12.0, 18.0],
            tau=0.05,
            paper_ready=True,
        ),
        "grid_ou_phi_asym": dict(
            mode="grid",
            recovery="asym",
            env_model="ou",
            fluct_target="Phi",
            env_sigma_list=[0.0, 0.10, 0.20, 0.30],
            k_list=[0.0, 0.5, 1.0, 1.5],  # multiplier on (k_minus,k_plus)
            tau=0.05,
            paper_ready=True,
        ),
    }


def run_preset(
    *,
    name: str,
    cfg: Config,
    eta_crit_list: List[float],
    outdir: Path,
    annotate: str = "none",
) -> Path:
    presets = preset_table()
    if name not in presets:
        raise ValueError(f"Unknown preset '{name}'. Use --list-presets to see available.")
    p = presets[name]
    mode = p["mode"]

    outdir.mkdir(parents=True, exist_ok=True)

    if mode == "overlay":
        outname = f"Rstruct_overlay_{name}.png"
        plot_overlay_three_stage(
            cfg=cfg,
            eta_crit_list=eta_crit_list,
            recovery=p["recovery"],
            env_model=p["env_model"],
            fluct_target=p["fluct_target"],
            env_sigma=float(p["env_sigma"]),
            tau=float(p["tau"]),
            outpath=outdir / outname,
            paper_ready=bool(p["paper_ready"]),
            show_textbox=True,
            textbox_mode="compact",
        )
        return outdir / outname

    if mode == "grid":
        outname = f"Rstruct_grid_{name}.png"
        plot_grid_recovery_vs_env(
            cfg=cfg,
            eta_crit_list=eta_crit_list,
            recovery=p["recovery"],
            env_model=p["env_model"],
            fluct_target=p["fluct_target"],
            env_sigma_list=list(p["env_sigma_list"]),
            k_list=list(p["k_list"]),
            tau=float(p["tau"]),
            outpath=outdir / outname,
            paper_ready=bool(p["paper_ready"]),
            annotate=annotate if annotate in ("none", "int", "2dp") else "none",
            remove_panel_legends=True,
            big_outer_labels=True,
        )
        return outdir / outname

    raise ValueError(f"Bad preset mode: {mode}")


# ============================================================
# CLI
# ============================================================

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser()

    # modes
    ap.add_argument("--overlay", action="store_true", help="Generate a 3-stage overlay plot.")
    ap.add_argument("--grid", action="store_true", help="Generate a paper-ready grid plot.")
    ap.add_argument("--preset", type=str, default=None, help="Run a named paper preset.")
    ap.add_argument("--run-all-presets", action="store_true", help="Generate one figure per preset (final test set).")
    ap.add_argument("--list-presets", action="store_true", help="List available presets and exit.")

    ap.add_argument("--outdir", type=str, default=".", help="Output directory for figures.")

    # core config overrides
    ap.add_argument("--N", type=int, default=Config.N)
    ap.add_argument("--T", type=float, default=Config.T)
    ap.add_argument("--dt", type=float, default=Config.dt)
    ap.add_argument("--eta0", type=float, default=Config.eta0)

    ap.add_argument("--lambda-min", type=float, default=Config.lambda_min)
    ap.add_argument("--lambda-max", type=float, default=Config.lambda_max)

    ap.add_argument("--gamma", type=float, default=Config.gamma)
    ap.add_argument("--Phi0", type=float, default=Config.Phi0)
    ap.add_argument("--alpha", type=float, default=Config.alpha)

    # recovery
    ap.add_argument("--recovery", type=str, choices=["none", "quad", "asym"], default="quad")
    ap.add_argument("--eta-star", type=float, default=Config.eta_star)
    ap.add_argument("--k", type=float, default=Config.k)
    ap.add_argument("--k-minus", type=float, default=Config.k_minus)
    ap.add_argument("--k-plus", type=float, default=Config.k_plus)

    # fluctuations
    ap.add_argument("--env-model", type=str, choices=["none", "white", "ou"], default="white")
    ap.add_argument("--fluct-target", type=str, choices=["lambda", "Phi"], default="lambda")
    ap.add_argument("--env-sigma", type=float, default=Config.env_sigma)
    ap.add_argument("--tau", type=float, default=Config.tau)

    # thresholds
    ap.add_argument(
        "--eta-crit-list",
        type=str,
        default=None,
        help="Comma-separated eta_crit list. Default: 10 points in [0.50,0.95].",
    )

    # grid sweeps
    ap.add_argument("--env-sigma-list", type=str, default="0,0.1,0.2,0.3",
                    help="(grid) Comma-separated env_sigma list.")
    ap.add_argument("--k-list", type=str, default="0,6,12,18",
                    help="(grid) For quad: k values. For asym: multipliers for (k_minus,k_plus).")

    ap.add_argument("--grid-annotate", type=str, default="none", choices=["none", "int", "2dp"],
                    help="Annotate points inside grid panels (none|int|2dp).")
    ap.add_argument("--grid-draft", action="store_true",
                    help="Use larger draft styling for grid (not compact paper-ready).")

    # RNG
    ap.add_argument("--seed", type=int, default=Config.seed)
    ap.add_argument("--no-crn", action="store_true", help="Disable common random numbers across conditions.")

    # overlay textbox controls
    ap.add_argument("--no-textbox", action="store_true",
                    help="Disable the parameter textbox on overlay plots.")
    ap.add_argument("--textbox", type=str, default="compact", choices=["compact", "full"],
                    help="Overlay textbox content. Default: compact.")

    return ap


def main() -> int:
    ap = build_parser()
    args = ap.parse_args()

    if args.list_presets:
        names = sorted(preset_table().keys())
        print("Available presets:")
        for n in names:
            print("  -", n)
        return 0

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    eta_crit_list = _parse_csv_floats(args.eta_crit_list) if args.eta_crit_list else _eta_crit_default()

    cfg = Config(
        N=int(args.N),
        T=float(args.T),
        dt=float(args.dt),
        eta0=float(args.eta0),
        lambda_min=float(args.lambda_min),
        lambda_max=float(args.lambda_max),
        gamma=float(args.gamma),
        Phi0=float(args.Phi0),
        alpha=float(args.alpha),
        eta_star=float(args.eta_star),
        k=float(args.k),
        k_minus=float(args.k_minus),
        k_plus=float(args.k_plus),
        env_sigma=float(args.env_sigma),
        tau=float(args.tau),
        seed=int(args.seed),
        crn=not bool(args.no_crn),
    )

    # If preset is specified, run it and exit (unless run-all-presets is also specified)
    if args.preset and not args.run_all_presets:
        path = run_preset(
            name=args.preset,
            cfg=cfg,
            eta_crit_list=eta_crit_list,
            outdir=outdir,
            annotate=str(args.grid_annotate),
        )
        print(f"[write] {path}")
        return 0

    # run all presets = final test set
    if args.run_all_presets:
        for name in sorted(preset_table().keys()):
            path = run_preset(
                name=name,
                cfg=cfg,
                eta_crit_list=eta_crit_list,
                outdir=outdir,
                annotate=str(args.grid_annotate),
            )
            print(f"[write] {path}")
        return 0

    # If neither overlay nor grid chosen, default to overlay
    do_overlay = bool(args.overlay) or (not args.overlay and not args.grid and not args.preset)
    do_grid = bool(args.grid)

    recovery: RecoveryType = str(args.recovery)  # type: ignore
    env_model: EnvModel = str(args.env_model)    # type: ignore
    fluct_target: FluctTarget = str(args.fluct_target)  # type: ignore

    if do_overlay:
        fname = f"Rstruct_overlay_rec{recovery}_env{env_model}_tgt{fluct_target}_sig{cfg.env_sigma:g}.png"
        plot_overlay_three_stage(
            cfg=cfg,
            eta_crit_list=eta_crit_list,
            recovery=recovery,
            env_model=env_model,
            fluct_target=fluct_target,
            env_sigma=cfg.env_sigma,
            tau=cfg.tau,
            outpath=outdir / fname,
            paper_ready=True,
            show_textbox=(not bool(args.no_textbox)),
            textbox_mode=str(args.textbox),
        )
        print(f"[write] {outdir / fname}")

    if do_grid:
        env_sigma_list = _parse_csv_floats(args.env_sigma_list)
        k_list = _parse_csv_floats(args.k_list)
        fname = f"Rstruct_grid_rec{recovery}_env{env_model}_tgt{fluct_target}.png"
        plot_grid_recovery_vs_env(
            cfg=cfg,
            eta_crit_list=eta_crit_list,
            recovery=recovery,
            env_model=env_model,
            fluct_target=fluct_target,
            env_sigma_list=env_sigma_list,
            k_list=k_list,
            tau=cfg.tau,
            outpath=outdir / fname,
            paper_ready=(not args.grid_draft),
            annotate=str(args.grid_annotate),  # type: ignore
            remove_panel_legends=True,
            big_outer_labels=True,
        )
        print(f"[write] {outdir / fname}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
