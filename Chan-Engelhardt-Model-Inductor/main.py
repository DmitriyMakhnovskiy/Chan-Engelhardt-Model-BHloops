"""
Chan-Engelhardt magnetic hysteresis model coupled to a single-loop
magnetic circuit with an air gap.

The simulation starts at t = 0 from the fully demagnetized state:

    H(0) = 0
    B(0) = 0
    i(0) = 0

All numerical parameters are stored only in Parameters.txt.
The voltage function V(t) is defined in V.py.
The hysteresis model is implemented in BH.py and is not modified.

Outputs:
    Current.csv
        time_s, I_A

    BHloop.csv
        H_A_per_m, B_T

The program also plots:
    i(t)
    the complete B(H) trajectory

Dependencies:
    Python 3
    NumPy
    Matplotlib

Dmitriy Makhnovskiy, 04.09.2026
"""

import numpy as np
import matplotlib.pyplot as plt
from BH import BH_increasing, BH_decreasing
from Parameters import parameters as par
from V import V


# ============================================================
# CONSTANTS
# ============================================================

pi = 3.1415926535897932384626433832795
mu0 = 4.0 * pi * 1.0e-7


# ============================================================
# TIME GRID
# ============================================================

n_steps = int(round(par["simulation_time"] / par["dt"]))

if not np.isclose(
    n_steps * par["dt"],
    par["simulation_time"]
):
    raise ValueError(
        "simulation_time must be an integer multiple of dt"
    )

t = np.arange(n_steps + 1, dtype=float) * par["dt"]


# ============================================================
# INITIAL FULLY DEMAGNETIZED STATE
# ============================================================

Hlast = 0.0
Blast = 0.0
Ilast = 0.0

OnInitMag = True
IsLatched = False

Bd = 0.0
Bd_turn = 0.0

H0 = 0.0
scale = 0.0
offset = 0.0

previous_direction = None


# ============================================================
# ARRAYS FOR RESULTS
# ============================================================

n_points = len(t)

Ivalues = np.zeros(n_points)
Hvalues = np.zeros(n_points)
Bvalues = np.zeros(n_points)

Ivalues[0] = 0.0
Hvalues[0] = 0.0
Bvalues[0] = 0.0


# ============================================================
# TRIAL HYSTERESIS EVALUATION
# ============================================================

def evaluate_hysteresis_trial(Htrial):

    if abs(Htrial - Hlast) <= par["H_tol"]:

        trial_state = (
            OnInitMag,
            IsLatched,
            Bd,
            Bd_turn,
            H0,
            scale,
            offset
        )

        return Blast, trial_state, previous_direction

    if Htrial > Hlast:

        if previous_direction is None:
            was_inc = True
        else:
            was_inc = (previous_direction == "inc")

        result = BH_increasing(
            Htrial,
            was_inc,
            Hlast,
            Blast,
            OnInitMag,
            IsLatched,
            Bd,
            Bd_turn,
            H0,
            scale,
            offset
        )

        Btrial = result[0]
        trial_state = result[1:]

        return Btrial, trial_state, "inc"

    if previous_direction is None:
        was_dec = True
    else:
        was_dec = (previous_direction == "dec")

    result = BH_decreasing(
        Htrial,
        was_dec,
        Hlast,
        Blast,
        OnInitMag,
        IsLatched,
        Bd,
        Bd_turn,
        H0,
        scale,
        offset
    )

    Btrial = result[0]
    trial_state = result[1:]

    return Btrial, trial_state, "dec"


# ============================================================
# RESIDUAL
# ============================================================

def residual(Htrial, Vnext):

    Btrial, trial_state, trial_direction = evaluate_hysteresis_trial(Htrial)

    Itrial = (
        par["l_Fe"] * Htrial
        + (par["g"] / mu0) * Btrial
    ) / par["N"]

    F = (
        par["R"] * Itrial
        + par["N"] * par["A"]
        * (Btrial - Blast) / par["dt"]
        - Vnext
    )

    return F, Btrial, Itrial, trial_state, trial_direction


# ============================================================
# SOLVE ONE TIME STEP
# ============================================================

def solve_time_step(Vnext):

    F0, B0, I0, state0, direction0 = residual(Hlast, Vnext)

    if not np.isfinite(F0):
        raise RuntimeError(
            "Non-finite residual at Hlast."
        )

    if abs(F0) <= par["residual_tol"]:
        return Hlast, B0, I0, state0, direction0

    # First try the physically expected direction.
    # If it fails to bracket the root, try the opposite direction.
    drive = Vnext - par["R"] * Ilast

    if drive > 0.0:
        primary_direction = 1.0
    elif drive < 0.0:
        primary_direction = -1.0
    else:
        if previous_direction == "dec":
            primary_direction = -1.0
        else:
            primary_direction = 1.0

    Ha = Hlast
    Fa = F0

    bracket_found = False

    for search_direction in (
        primary_direction,
        -primary_direction
    ):

        step = par["initial_H_step"]

        for _ in range(int(par["max_bracket_iter"])):

            Hb = Hlast + search_direction * step
            Fb, _, _, _, _ = residual(Hb, Vnext)

            if not np.isfinite(Fb):
                raise RuntimeError(
                    "Non-finite residual while searching for H bracket."
                )

            if Fa * Fb <= 0.0:
                bracket_found = True
                break

            step *= 2.0

        if bracket_found:
            break

    if not bracket_found:
        raise RuntimeError(
            "Unable to bracket H root in either direction. "
            "Reduce dt or check circuit/material parameters."
        )

    if Ha > Hb:
        Ha, Hb = Hb, Ha
        Fa, Fb = Fb, Fa

    best = None
    best_abs_F = np.inf

    for _ in range(int(par["max_bisection_iter"])):

        Hm = 0.5 * (Ha + Hb)

        (
            Fm,
            Bm,
            Im,
            state_m,
            direction_m
        ) = residual(Hm, Vnext)

        if not np.isfinite(Fm):
            raise RuntimeError(
                "Non-finite residual during bisection."
            )

        if abs(Fm) < best_abs_F:
            best_abs_F = abs(Fm)
            best = (Hm, Bm, Im, state_m, direction_m)

        if (
            abs(Fm) <= par["residual_tol"]
            or abs(Hb - Ha) <= par["H_tol"]
        ):
            return best

        if Fa * Fm <= 0.0:
            Hb = Hm
            Fb = Fm
        else:
            Ha = Hm
            Fa = Fm

    if best is None:
        raise RuntimeError("Bisection solver failed unexpectedly.")

    raise RuntimeError(
        "Bisection did not converge within max_bisection_iter. "
        f"Best |residual| = {best_abs_F:.6e}."
    )


# ============================================================
# TIME INTEGRATION
# ============================================================

for n in range(n_points - 1):

    Vnext = V(t[n + 1])

    (
        Hnext,
        Bnext,
        Inext,
        state_next,
        direction_next
    ) = solve_time_step(Vnext)

    (
        OnInitMag,
        IsLatched,
        Bd,
        Bd_turn,
        H0,
        scale,
        offset
    ) = state_next

    if abs(Hnext - Hlast) > par["H_tol"]:
        previous_direction = direction_next

    Ivalues[n + 1] = Inext
    Hvalues[n + 1] = Hnext
    Bvalues[n + 1] = Bnext

    Hlast = Hnext
    Blast = Bnext
    Ilast = Inext


# ============================================================
# SAVE CURRENT VERSUS TIME
# ============================================================

np.savetxt(
    "Current.csv",
    np.column_stack((t, Ivalues)),
    delimiter=",",
    header="time_s,I_A",
    comments=""
)


# ============================================================
# SAVE COMPLETE B-H TRAJECTORY
# ============================================================

np.savetxt(
    "BHloop.csv",
    np.column_stack((Hvalues, Bvalues)),
    delimiter=",",
    header="H_A_per_m,B_T",
    comments=""
)


# ============================================================
# PLOT CURRENT VERSUS TIME
# ============================================================

plt.figure(figsize=(8, 5))
plt.plot(t, Ivalues, linewidth=1.5)
plt.xlabel("Time (s)")
plt.ylabel("Current (A)")
plt.title("Coil Current")
plt.grid(True)
plt.tight_layout()


# ============================================================
# PLOT COMPLETE B-H TRAJECTORY
# ============================================================

plt.figure(figsize=(7, 6))
plt.plot(Hvalues, Bvalues, linewidth=1.5)
plt.xlabel("H (A/m)")
plt.ylabel("B (T)")
plt.title("B-H Hysteresis Trajectory")
plt.grid(True)
plt.axhline(0.0, linewidth=0.8)
plt.axvline(0.0, linewidth=0.8)
plt.tight_layout()

plt.show()
