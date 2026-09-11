"""
External voltage source V(t).

Define the required voltage waveform in the function V(t).

All numerical parameters used by the voltage function should be stored
in Parameters.txt and accessed through par["parameter_name"].
"""

import numpy as np
from Parameters import parameters as par


def V(t):
    return (
        par["V_offset"]
        + par["V_amplitude"]
        * np.sin(2.0 * np.pi * par["V_frequency"] * t)
    )
