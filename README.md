# CCST-structural-stability-recovery-env
Stability Geometry of Branch Structures in Open Quantum Systems (Extended): Recovery Geometry and Environmental Coupling

% ============================================================
% Appendix A (README) — ccst_env_recovery_sweep_vnext.py
% Copy/paste-ready LaTeX block (FINAL)
%  - Matches the user's current script exactly.
%  - Includes the small “alignment fixes” discussed:
%      (i) env-model default clarified (without presets),
%     (ii) recovery=none makes grid row-sweep redundant,
%    (iii) asym grid uses multipliers; s=0 disables recovery.
% ============================================================

\section*{Appendix A. README: \texttt{ccst\_env\_recovery\_sweep\_vnext.py}}

\subsection*{A.1 Purpose and scope (honesty statement)}

\texttt{ccst\_env\_recovery\_sweep\_vnext.py} is a \textbf{toy-model figure generator} for a reduced,
bounded proxy state $\eta(t)\in[0,1]$ that represents \emph{structural coherence / stability} in the
manuscript’s stability-geometry discussion.  Its purpose is to produce \textbf{paper-ready, reproducible}
visualizations that cleanly separate and illustrate:

\begin{enumerate}
\item \textbf{Baseline environment coupling} (the deterministic drift backbone),
\item \textbf{Recovery geometry} (restoring drift; symmetric \texttt{quad} or direction-asymmetric \texttt{asym}),
\item \textbf{Environmental fluctuations} (fast \texttt{white} noise or correlated \texttt{ou} noise),
      acting on either:
      \begin{itemize}
      \item \textbf{target=\texttt{lambda}} (effective coupling rate channel), or
      \item \textbf{target=\texttt{Phi}} (environmental drive channel).
      \end{itemize}
\end{enumerate}

\paragraph{Honesty note (what this is not).}
This script is \emph{not} a microscopic branching model.  It does not define probabilities, Born weights,
or a selection rule.  It is intentionally a \textbf{visual regime toy}: it isolates a deterministic
drift backbone, then adds (i) recovery geometry and (ii) environment variability as fluctuations acting
on a specified target channel.  Its value is interpretability and auditability, not microphysical fidelity.

\subsection*{A.2 Model definition}

We simulate an ensemble of $N$ trajectories of a bounded proxy state $\eta_b(t)\in[0,1]$ for branch label $b$.

\paragraph{Baseline deterministic backbone (environment-only).}
Each trajectory draws a heterogeneous coupling parameter $\lambda_b$ once at initialization:
\[
\lambda_b \sim U(\lambda_{\min},\lambda_{\max}).
\]
The environment-only deterministic dynamics are:
\begin{equation}
\frac{d\eta_b}{dt} \;=\; -\lambda_b\,\eta_b \;-\; \gamma\,\Phi(\eta_b),
\qquad
\Phi(\eta) \;=\; \Phi_0 + \alpha(1-\eta).
\end{equation}

\paragraph{Recovery geometry (optional).}
Recovery adds a restoring drift toward a target level $\eta^\ast$:
\begin{equation}
\frac{d\eta_b}{dt}
\;=\;
-\lambda_b\,\eta_b
-\gamma\,\Phi(\eta_b)
+\;k_{\mathrm{side}}(\eta_b)\,\big(\eta^\ast-\eta_b\big).
\end{equation}
The script supports:
\begin{itemize}
\item \texttt{recovery=none}: $k_{\mathrm{side}}(\eta)=0$.
\item \texttt{recovery=quad}: $k_{\mathrm{side}}(\eta)=k$ (symmetric recovery strength).
\item \texttt{recovery=asym}:
\[
k_{\mathrm{side}}(\eta)=
\begin{cases}
k_- & \eta < \eta^\ast,\\
k_+ & \eta \ge \eta^\ast.
\end{cases}
\]
This implements direction-asymmetric recovery while keeping the same basin center $\eta^\ast$.
\end{itemize}

\subsection*{A.3 Environmental fluctuations (optional)}

Fluctuations can act on \textbf{$\lambda$} or \textbf{$\Phi$}.

\begin{itemize}
\item \textbf{target=\texttt{lambda}}:
\[
\lambda_b \rightarrow \lambda_b + \delta\lambda_b(t)
\quad\Rightarrow\quad
\text{extra drift } -\delta\lambda_b(t)\,\eta_b.
\]
\item \textbf{target=\texttt{Phi}}:
\[
\Phi \rightarrow \Phi + \delta\Phi_b(t)
\quad\Rightarrow\quad
\text{extra drift } -\gamma\,\delta\Phi_b(t).
\]
\end{itemize}

Two fluctuation models are supported:

\paragraph{(i) \texttt{white} noise (Euler--Maruyama increment).}
At each time step, an Euler--Maruyama (EM) increment is applied:
\begin{itemize}
\item target=\texttt{lambda}: \;\; $\eta \leftarrow \eta - (\sigma_{\mathrm{env}}\,\eta)\sqrt{\Delta t}\,Z$
\item target=\texttt{Phi}: \;\;\;\;\;\;\; $\eta \leftarrow \eta - (\gamma\,\sigma_{\mathrm{env}})\sqrt{\Delta t}\,Z$
\end{itemize}
where $Z\sim\mathcal{N}(0,1)$.  In this mode, $\sigma_{\mathrm{env}}$ is a diffusion amplitude (per
$\sqrt{\text{time}}$).  Equivalently, this corresponds to an It\^{o} SDE containing the stated diffusion term,
discretized by Euler--Maruyama \emph{together with} the deterministic drift terms above.

\paragraph{(ii) \texttt{ou} noise (correlated fluctuations).}
An OU state $X(t)$ is evolved with correlation time $\tau$ and stationary standard deviation $\sigma_{\mathrm{env}}$:
\begin{equation}
dX = -\frac{1}{\tau}X\,dt + \sqrt{\frac{2}{\tau}}\sigma_{\mathrm{env}}\,dW.
\end{equation}
The script uses an exact OU discretization for the OU update:
\[
X_{t+\Delta t} = aX_t + \sigma_{\mathrm{env}}\sqrt{1-a^2}\,Z,\quad a=e^{-\Delta t/\tau},
\]
so that the OU stationary variance is preserved (i.e., $\mathrm{Std}(X)=\sigma_{\mathrm{env}}$).
After updating $X$, it is inserted into the drift channel:
\begin{itemize}
\item target=\texttt{lambda}: \;\; add drift $-\,(X\,\eta)$
\item target=\texttt{Phi}: \;\;\;\;\;\;\; add drift $-\gamma X$
\end{itemize}

\paragraph{Bounding.}
After each update, $\eta$ is clipped to $[0,1]$. This enforces the interpretation of $\eta$ as a bounded proxy.

\subsection*{A.4 Diagnostic: structural survival fraction $R_{\mathrm{struct}}$}

For a threshold $\eta_{\mathrm{crit}}$, a trajectory survives on $[0,T]$ iff
\[
\min_{t\in[0,T]} \eta(t) \ge \eta_{\mathrm{crit}}.
\]
The structural survival fraction is:
\begin{equation}
R_{\mathrm{struct}}(\eta_{\mathrm{crit}})
=
\frac{1}{N}\sum_{b=1}^{N}
\mathbf{1}\!\left[\min_{t\in[0,T]}\eta_b(t)\ge \eta_{\mathrm{crit}}\right].
\end{equation}
All main figures plot $R_{\mathrm{struct}}$ (typically in \%) as a function of $\eta_{\mathrm{crit}}$.

\subsection*{A.5 Figure types produced by the script}

\paragraph{(A) Overlay (three-stage) plot.}
A single plot with three curves:
\begin{enumerate}
\item env-only: no recovery, no fluctuations,
\item + recovery: recovery enabled, no fluctuations,
\item + recovery + fluctuations: recovery enabled and fluctuations enabled.
\end{enumerate}
This isolates the incremental effect of (i) recovery geometry and (ii) environmental variability.

\paragraph{(B) Grid plot (paper-ready).}
A panel grid with:
\[
\text{columns}=\sigma_{\mathrm{env}},\qquad
\text{rows}=\text{recovery strength}.
\]
Panels share identical y-axis limits (shared scale) so that \emph{small differences remain visually comparable}.
The grid is designed to avoid panel-level clutter (no per-panel legends) and instead uses outer labels.

\subsection*{A.6 Reproducibility: RNG, seeds, and CRN}

The script supports stable runs via \verb|--seed|.

\paragraph{CRN (common random numbers).}
By default \verb|crn=True|, meaning \textbf{common random numbers} are used across conditions to stabilize
comparisons (reducing Monte Carlo variance in differences). Disable via:
\begin{verbatim}
--no-crn
\end{verbatim}
Disabling CRN is useful to test robustness: qualitative patterns should persist without variance reduction.

\subsection*{A.7 Command-line usage (complete)}

\paragraph{Modes.}
\begin{itemize}
\item \verb|--overlay| : generate a 3-stage overlay plot.
\item \verb|--grid| : generate a paper-ready grid plot.
\item If neither is specified, the default behavior is to generate an overlay.
\end{itemize}

\paragraph{Presets.}
\begin{itemize}
\item \verb|--preset <name>| : run one named preset configuration (overlay or grid).
\item \verb|--run-all-presets| : generate one figure per preset (``one of each kind'' final test set).
\item \verb|--list-presets| : list available preset names.
\end{itemize}

\paragraph{Core parameters.}
\begin{itemize}
\item \verb|--N| : ensemble size ($N$ trajectories).
\item \verb|--T| : time horizon $T$.
\item \verb|--dt| : step size $\Delta t$.
\item \verb|--eta0| : initial condition $\eta(0)$.
\end{itemize}

\paragraph{Baseline coupling and drive.}
\begin{itemize}
\item \verb|--lambda-min|, \verb|--lambda-max| : define the sampling range for $\lambda_b$.
\item \verb|--gamma| : strength of the $\Phi$ channel.
\item \verb|--Phi0| : base drive level.
\item \verb|--alpha| : state dependence of $\Phi(\eta)=\Phi_0+\alpha(1-\eta)$.
\end{itemize}

\paragraph{Recovery geometry.}
\begin{itemize}
\item \verb|--recovery {none,quad,asym}|.
\item \verb|--eta-star| sets $\eta^\ast$.
\item For \texttt{quad}: \verb|--k| sets the restoring strength $k$.
\item For \texttt{asym}: \verb|--k-minus| and \verb|--k-plus| set direction-dependent strengths $(k_-,k_+)$.
\end{itemize}

\paragraph{Fluctuations.}
\begin{itemize}
\item \verb|--env-model {none,white,ou}| (default: \texttt{white} when not using presets).
\item \verb|--fluct-target {lambda,Phi}|.
\item \verb|--env-sigma| : fluctuation amplitude. Fluctuations are effectively ``off'' if \verb|--env-sigma 0|.
\item \verb|--tau| : OU correlation time (only used if \texttt{ou}).
\end{itemize}

\paragraph{Threshold grid.}
\begin{itemize}
\item \verb|--eta-crit-list| : comma-separated list. If omitted, the script uses
10 evenly spaced points in $[0.50,0.95]$.
\end{itemize}

\paragraph{Grid sweep lists.}
\begin{itemize}
\item \verb|--env-sigma-list| : comma-separated values for the grid columns.
\item \verb|--k-list| : row values:
  \begin{itemize}
  \item for \texttt{quad}: literal $k$ values,
  \item for \texttt{asym}: multiplicative scale factors $s$ applied as $(k_-,k_+)\leftarrow s\,(k_-,k_+)$.
        In particular, $s=0$ turns recovery off.
  \end{itemize}
\item If \verb|--recovery none|, the grid row sweep over \verb|--k-list| is redundant (rows become identical up to Monte Carlo noise).
\item \verb|--grid-annotate {none,int,2dp}| : optional annotation inside panels (normally \texttt{none} for paper).
\item \verb|--grid-draft| : use larger draft styling for the grid (helpful for debugging).
\end{itemize}

\subsection*{A.8 Presets included in \texttt{ccst\_env\_recovery\_sweep\_vnext.py}}

The script defines presets in \verb|preset_table()|. At the time of writing they include:

\begin{itemize}
\item \verb|overlay_white_lambda_quad|
\item \verb|overlay_white_lambda_asym|
\item \verb|overlay_ou_phi_quad|
\item \verb|overlay_ou_phi_asym|
\item \verb|grid_white_lambda_quad|
\item \verb|grid_white_lambda_asym|
\item \verb|grid_ou_phi_quad|
\item \verb|grid_ou_phi_asym|
\end{itemize}

\paragraph{Example.}
\begin{verbatim}
python3 ccst_env_recovery_sweep_vnext.py --preset grid_ou_phi_asym --outdir figs
python3 ccst_env_recovery_sweep_vnext.py --run-all-presets --outdir figs
\end{verbatim}

\subsection*{A.9 Practical parameter intuition (how to explore)}

\paragraph{Common ways to decrease survival.}
\begin{itemize}
\item increase \verb|--lambda-max| (stronger coupling),
\item increase \verb|--gamma| or \verb|--Phi0| (stronger $\Phi$ channel),
\item increase \verb|--env-sigma| (stronger fluctuations),
\item increase \verb|--T| (more time to degrade),
\item decrease \verb|--eta0| (worse initial condition).
\end{itemize}

\paragraph{Common ways to increase survival.}
\begin{itemize}
\item increase recovery strength (\verb|--k| or \verb|--k-minus|/\verb|--k-plus|),
\item choose $\eta^\ast$ within the high-coherence region to represent a stable basin center.
\end{itemize}

\paragraph{Asymmetry guidance.}
If the intended narrative is ``fast recovery from degraded states, slower erosion otherwise'',
a typical choice is $k_- > k_+$.

\subsection*{A.10 Troubleshooting}

\paragraph{Curves look too similar.}
Try increasing \verb|--env-sigma|, increasing \verb|--lambda-max|, increasing $T$,
widening the recovery sweep range in \verb|--k-list|, or increasing $N$.
If the goal is clean comparisons, keep CRN on; if the goal is robustness,
turn CRN off and confirm qualitative stability.

\paragraph{Numerical stability.}
If artifacts appear, reduce \verb|--dt| and verify that results are stable under timestep refinement.

% ============================================================
% End Appendix A README (FINAL)
% ============================================================
