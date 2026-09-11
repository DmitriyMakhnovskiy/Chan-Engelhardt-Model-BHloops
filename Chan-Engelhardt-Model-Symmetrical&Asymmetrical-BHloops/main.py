
"""
Chan-Engelhardt magnetic hysteresis model.

This program calculates B-H hysteresis loops for a prescribed magnetic-field
history H(t). The model supports major and minor hysteresis loops, including
asymmetrical nested loops and arbitrary field reversals.

Project structure:
    main.py
        Main simulation script. Generates the H(t) excitation history,
        calculates the corresponding B(H) trajectory, plots the hysteresis
        loop, and saves the result to BHloop.csv.

    BH.py
        Core hysteresis model. Contains the major-loop functions, initial
        magnetization curve, reversal logic, and the BH_increasing() and
        BH_decreasing() routines used to construct history-dependent
        hysteresis branches.

    Hscan.py
        Defines test excitation functions H(t), including damped and biased
        waveforms for generating symmetrical and asymmetrical minor loops.

    Parameters.py
        Reads magnetic material parameters from Parameters.txt and provides
        them to the other program modules.

    Parameters.txt
        Contains the magnetic material parameters:
            Bs - saturation magnetic flux density [T]
            Br - remanent magnetic flux density [T]
            Hc - coercive magnetic field strength [A/m]

    BHloop.csv
        Output file containing the calculated hysteresis trajectory in two
        columns:
            H_A_per_m, B_T

Dependencies:
    Python 3
    NumPy
    Matplotlib

Dmitriy Makhnovskiy, 04.09.2026

"""

import numpy as np
import matplotlib.pyplot as plt
from BH import BH_increasing, BH_decreasing
from Hscan import Hscan1, Hscan2


# ============================================================
# GENERATE H(t)
# ============================================================

t_max = 0.1  # time, s
t = np.linspace(0, t_max, 1000)

# Choose a function Hscan(1,2,...) from Hscan.py
Hvalues = Hscan1(t)

# ============================================================
# CHECK INITIAL CONDITION
# ============================================================

if not np.isclose(Hvalues[0], 0.0):
    print('Error: simulation must start with H = 0')
    raise SystemExit


# ============================================================
# INITIAL STATE
# ============================================================

Hlast = 0.0
Blast = 0.0

OnInitMag = True
IsLatched = False

Bd = 0.0
Bd_turn = 0.0

H0 = 0.0
scale = 0.0
offset = 0.0


# ============================================================
# ARRAY FOR B VALUES
# ============================================================

Bvalues = np.zeros(len(Hvalues))

# Initial point:
# H(0) = 0
# B(0) = 0
Bvalues[0] = 0.0


# ============================================================
# PREVIOUS DIRECTION
#
# "inc"  - H was increasing
# "dec"  - H was decreasing
# None   - no direction has yet been established
# ============================================================

previous_direction = None


# ============================================================
# CALCULATE B(H)
# ============================================================

for i, H in enumerate(Hvalues[1:], start=1):

    # ========================================================
    # H HAS NOT CHANGED
    # ========================================================

    if np.isclose(H, Hlast):

        # B and all state variables remain unchanged.
        B = Blast

        # Important:
        # previous_direction is NOT changed.


    # ========================================================
    # CURRENT H IS INCREASING
    # ========================================================

    elif H > Hlast:

        # On the first movement there was no reversal.
        if previous_direction is None:
            was_inc = True
        else:
            was_inc = (previous_direction == "inc")

        (
            B,
            OnInitMag,
            IsLatched,
            Bd,
            Bd_turn,
            H0,
            scale,
            offset
        ) = BH_increasing(
            H,
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

        previous_direction = "inc"


    # ========================================================
    # CURRENT H IS DECREASING
    # ========================================================

    else:

        # On the first movement there was no reversal.
        if previous_direction is None:
            was_dec = True
        else:
            was_dec = (previous_direction == "dec")

        (
            B,
            OnInitMag,
            IsLatched,
            Bd,
            Bd_turn,
            H0,
            scale,
            offset
        ) = BH_decreasing(
            H,
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

        previous_direction = "dec"


    # ========================================================
    # SAVE CURRENT B VALUE
    # ========================================================

    Bvalues[i] = B


    # ========================================================
    # CURRENT POINT BECOMES THE PREVIOUS POINT
    # ========================================================

    Hlast = H
    Blast = B


# ============================================================
# CREATE B-H ARRAY
# ============================================================

BHloop = np.column_stack((Hvalues, Bvalues))


# ============================================================
# DRAW B(H)
# ============================================================

plt.figure(figsize=(7, 6))

plt.plot(
    BHloop[:, 0],
    BHloop[:, 1],
    linewidth=1.5
)

plt.xlabel('H (A/m)')
plt.ylabel('B (T)')
plt.title('B-H Hysteresis Loop')

plt.grid(True)

plt.axhline(0.0, linewidth=0.8)
plt.axvline(0.0, linewidth=0.8)

plt.tight_layout()
plt.show()

# ============================================================
# SAVE B-H LOOP TO CSV FILE
# ============================================================

np.savetxt(
    "BHloop.csv",
    BHloop,
    delimiter=",",
    header="H_A_per_m,B_T",
    comments=""
)