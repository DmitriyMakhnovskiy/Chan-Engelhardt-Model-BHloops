# Chan–Engelhardt Model for B–H Loops and Nonlinear Inductors

This repository contains Python and LTspice implementations of analytical magnetic-hysteresis models based on the **Chan model** and the **Chan–Engelhardt extension** for asymmetric minor hysteresis loops.

The project is intended for both educational and engineering use. It includes:

- calculation of major and minor B–H hysteresis loops;
- modelling of magnetic memory and reversal points;
- the Engelhardt correction for non-physical asymmetric minor loops produced by the original Chan approach;
- transient simulation of an inductor with a ferromagnetic core and optional air gap;
- implicit time integration using the **Backward Euler** method;
- nonlinear root bracketing and bisection;
- comparison of the Python transient solver with **LTspice**;
- supporting reports, papers, patent material, and introductory presentations on electrodynamics and magnetic circuits.

---

## Repository structure

```text
Chan-Engelhardt-Model-BHloops/
│
├── Chan-Model-Symmetrical-BHloops/
│   ├── Chan_model.py
│   ├── Gapped_BH-loop_data.csv
│   ├── Ungapped_BH-loop_data.csv
│   └── simulation_log.log
│
├── Chan-Engelhardt-Model-Symmetrical&Asymmetrical-BHloops/
│   ├── main.py
│   ├── BH.py
│   ├── Hscan.py
│   ├── Parameters.py
│   ├── Parameters.txt
│   ├── BHloop.csv
│   ├── Asymmetric Minor Hysteresis Loop Model and Circuit Simulator Including the Same_EN.pdf
│   ├── Asymmetric Minor Hysteresis Loop Model and Circuit Simulator Including the Same_RUS.pdf
│   ├── ASYMMETRIC MINOR HYSTERESIS LOOPS_US7502723.pdf
│   └── Jiles-Atherton model.pdf
│
├── Chan-Engelhardt-Model-Inductor/
│   ├── main.py
│   ├── BH.py
│   ├── V.py
│   ├── Parameters.py
│   ├── Parameters.txt
│   ├── Current.csv
│   ├── BHloop.csv
│   ├── LTspice_Inductor with a magnetic core/
│   ├── Asymmetric Minor Hysteresis Loop Model and Circuit Simulator Including the Same_EN.pdf
│   ├── Asymmetric Minor Hysteresis Loop Model and Circuit Simulator Including the Same_RUS.pdf
│   ├── ASYMMETRIC MINOR HYSTERESIS LOOPS_US7502723.pdf
│   └── Jiles-Atherton model.pdf
│
├── LTspice_Non-linear_Transformer/
│
├── A case study on magnetic cores_All.pdf
├── ASYMMETRIC MINOR HYSTERESIS LOOPS_US7502723.pdf
├── Nonlinear_transformer_model_for_circuit_simulation.pdf
├── Introduction to electrodynamics.ppsx
├── Magnetic circuits.ppsx
└── README.md
```

---

## 1. Chan Model: Symmetrical B–H Loops

Folder:

```text
Chan-Model-Symmetrical-BHloops/
```

This folder contains an implementation of the original analytical hysteresis model proposed by **Chan et al.**

The model provides a compact representation of the major hysteresis loop and symmetrical minor loops using the principal magnetic parameters:

- saturation flux density `Bs`;
- remanent flux density `Br`;
- coercive field strength `Hc`.

The main implementation is contained in:

```text
Chan_model.py
```

The folder also contains calculated data for gapped and ungapped magnetic circuits:

```text
Gapped_BH-loop_data.csv
Ungapped_BH-loop_data.csv
```

This project is useful as a reference implementation of the original Chan model before introducing the Engelhardt correction for asymmetric minor loops.

---

## 2. Chan–Engelhardt Model: Symmetrical and Asymmetrical B–H Loops

Folder:

```text
Chan-Engelhardt-Model-Symmetrical&Asymmetrical-BHloops/
```

This is the main Python implementation of the **Chan–Engelhardt hysteresis algorithm**.

The program reproduces both symmetrical and asymmetrical B–H trajectories while preserving the magnetic history required to describe field reversals and nested minor loops.

The implementation follows the logic of U.S. Patent **US7502723B1**.

### `BH.py`

Contains the magnetic material model, including:

- upper and lower branches of the major hysteresis loop;
- the initial magnetization curve;
- ordinary vertically shifted minor-loop branches;
- basis functions for the Engelhardt affine correction;
- limiting functions preventing corrected branches from crossing the preceding loop;
- numerical determination of `H0`;
- calculation of affine-map coefficients;
- `BH_increasing()` and `BH_decreasing()` state-dependent algorithms.

### `Hscan.py`

Defines prescribed magnetic-field histories `H(t)` used for testing:

- major loops;
- symmetrical minor loops;
- asymmetric minor loops;
- nested reversals;
- biased excitation histories.

### `main.py`

Controls the sequential scan through `H(t)`.

Its main functions are:

- preservation of the hysteresis state;
- determination of the direction of field variation;
- selection of the increasing or decreasing hysteresis branch;
- storage of the calculated B–H trajectory;
- plotting of the resulting hysteresis loop.

### `Parameters.py` and `Parameters.txt`

All numerical parameters are stored in the text file:

```text
Parameters.txt
```

and loaded by:

```text
Parameters.py
```

This keeps numerical parameters separate from the physical model.

### Output

The calculated hysteresis trajectory is written to:

```text
BHloop.csv
```

---

## 3. Chan–Engelhardt Model of a Nonlinear Inductor

Folder:

```text
Chan-Engelhardt-Model-Inductor/
```

This project couples the same Chan–Engelhardt hysteresis model to the electrical equations of an inductor with a ferromagnetic core and an optional air gap.

In this case the magnetic field `H(t)` is **not prescribed externally**.

Instead, it is determined self-consistently at every time step from:

- the applied voltage;
- coil resistance;
- magnetic-circuit geometry;
- the hysteretic material relation `B(H, state)`.

The transient problem is solved using an **implicit Backward Euler** time discretization.

At every time step, a nonlinear algebraic equation is solved for the new magnetic field.

### `BH.py`

Contains the same Chan–Engelhardt hysteresis model used in the standalone B–H-loop project.

### `V.py`

Defines the externally applied voltage waveform:

```text
V(t)
```

The voltage function can be modified independently of the nonlinear solver.

### `main.py`

Implements the transient nonlinear solver.

Its main functions are:

- construction of the time grid;
- preservation of the confirmed hysteresis state from the previous time step;
- trial evaluation of `B(H)` without altering that confirmed state;
- calculation of the nonlinear residual;
- automatic search for a root bracket;
- solution of the nonlinear equation by bisection;
- acceptance of the new hysteresis state only after convergence;
- calculation and storage of `I(t)`, `H(t)`, and `B(t)`.

A key numerical principle is that every trial point used during the nonlinear root search is evaluated from the **same confirmed hysteresis state**.

Trial states are not allowed to alter the magnetic history until the nonlinear equation has converged.

### `Parameters.py` and `Parameters.txt`

These files contain:

- magnetic material parameters;
- coil parameters;
- magnetic-core geometry;
- air-gap length;
- voltage-source parameters;
- simulation time;
- time step;
- nonlinear-solver tolerances;
- root-bracketing parameters;
- bisection parameters.

### Outputs

The calculated coil current is written to:

```text
Current.csv
```

The corresponding magnetic trajectory is written to:

```text
BHloop.csv
```

The program also plots:

```text
I(t)
B(H)
```

### LTspice comparison

The subfolder:

```text
LTspice_Inductor with a magnetic core/
```

contains LTspice simulations used for direct comparison with the Python transient solver.

The Python results were developed specifically to reproduce and verify the nonlinear-inductor behaviour obtained with the Chan–Engelhardt model in LTspice.

---

## 4. LTspice Nonlinear Transformer Example

Folder:

```text
LTspice_Non-linear_Transformer/
```

This folder contains an LTspice example of a nonlinear transformer using a magnetic-core model.

It is included as a reference for circuit-simulator implementations of nonlinear magnetic components and for comparison with the analytical Python models.

---

# Reports and Reference Material

## Main Technical Report

The most complete description of the project is provided in two language versions:

### English

```text
Asymmetric Minor Hysteresis Loop Model and Circuit Simulator Including the Same_EN.pdf
```

### Russian

```text
Asymmetric Minor Hysteresis Loop Model and Circuit Simulator Including the Same_RUS.pdf
```

The report covers:

- the physical and mathematical basis of the Chan model;
- the original Chan algorithm for symmetrical minor loops;
- the defect of the original treatment of some asymmetric minor loops;
- the Engelhardt correction introduced in US7502723B1;
- interpretation of the patent figures and flowcharts;
- construction of ordinary and affine-corrected minor branches;
- the complete hysteresis-state algorithm;
- the Python implementation;
- the state variables used to preserve magnetic memory;
- calculation of symmetrical and asymmetrical B–H trajectories;
- coupling of the hysteresis model to an electrical circuit;
- transient simulation of an inductor with a magnetic core;
- Backward Euler time discretization;
- nonlinear root bracketing;
- bisection;
- preservation of the magnetic state during nonlinear iterations;
- comparison of transient-current calculations with LTspice.

Copies of the reports are stored in the relevant Python project folders.

---

## Engelhardt Patent

```text
ASYMMETRIC MINOR HYSTERESIS LOOPS_US7502723.pdf
```

**M. Thomas Engelhardt**

*Asymmetric Minor Hysteresis Loop Model and Circuit Simulator Including the Same*

U.S. Patent **US7502723B1**, 2009.

The patent introduces the method used to correct problematic asymmetric minor loops in the original Chan model.

The corrected branch is constructed using an affine transformation of a basis function together with a limiting function that prevents the trajectory from producing the non-physical excursions outside the major hysteresis loop.

---

## Original Chan Paper

```text
Nonlinear_transformer_model_for_circuit_simulation.pdf
```

**J. H. Chan et al.**

*Nonlinear Transformer Model for Circuit Simulation*

IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems, Vol. 10, No. 4, 1991.

This paper provides the analytical hysteresis model on which the later Engelhardt extension is based.

---

## Jiles–Atherton Material

```text
Jiles-Atherton model.pdf
```

This material is included as background on the Jiles–Atherton approach to ferromagnetic hysteresis modelling.

The Jiles–Atherton model represents a different, more physically motivated approach to hysteresis and provides useful context for comparison with the compact analytical Chan model.

---

## Additional Magnetic-Core Report

```text
A case study on magnetic cores_All.pdf
```

Additional technical material related to magnetic cores and nonlinear magnetic behaviour.

---

## Introductory Presentations

### Introduction to Electrodynamics

```text
Introduction to electrodynamics.ppsx
```

Introductory material on electrodynamics used as theoretical background for the magnetic-circuit modelling work.

### Magnetic Circuits

```text
Magnetic circuits.ppsx
```

Engineering introduction to magnetic circuits and their analysis.

---

# Numerical Philosophy

The project deliberately separates the **magnetic material model** from the **electrical circuit solver**.

For a prescribed magnetic-field history, the standalone hysteresis project solves:

```text
H(t) -> B(t)
```

In the nonlinear-inductor model, the magnetic field is itself an unknown:

```text
V(t) -> H(t), B(t), I(t)
```

Because the hysteresis model contains magnetic memory, it should more accurately be regarded as:

```text
B = B(H, state)
```

rather than as an ordinary single-valued function:

```text
B = B(H)
```

This distinction is essential for the nonlinear inductor solver.

During the root search, every trial magnetic field generates a corresponding trial magnetic state, but the previously confirmed hysteresis state remains unchanged.

Only after the nonlinear equation has converged is the trial state accepted as the new physical state of the system.

---

# Requirements

The Python projects require:

```text
Python 3
NumPy
Matplotlib
```

Install the required packages with:

```bash
pip install numpy matplotlib
```

---

# Running the Standalone Hysteresis Model

Open the folder:

```text
Chan-Engelhardt-Model-Symmetrical&Asymmetrical-BHloops
```

Adjust the magnetic material parameters in:

```text
Parameters.txt
```

If required, modify the excitation history in:

```text
Hscan.py
```

Run:

```bash
python main.py
```

The program calculates and plots the B–H trajectory and writes the result to:

```text
BHloop.csv
```

---

# Running the Nonlinear-Inductor Model

Open the folder:

```text
Chan-Engelhardt-Model-Inductor
```

Adjust the model parameters in:

```text
Parameters.txt
```

These include:

- magnetic material parameters;
- coil parameters;
- core geometry;
- air-gap length;
- voltage parameters;
- time step;
- simulation duration;
- nonlinear-solver tolerances.

The voltage waveform is defined in:

```text
V.py
```

Run:

```bash
python main.py
```

The program generates:

```text
Current.csv
BHloop.csv
```

and plots:

```text
I(t)
B(H)
```

---

# References

1. **J. H. Chan et al.**,  
   “Nonlinear Transformer Model for Circuit Simulation,”  
   *IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems*,  
   Vol. 10, No. 4, 1991.

2. **M. Thomas Engelhardt**,  
   “Asymmetric Minor Hysteresis Loop Model and Circuit Simulator Including the Same,”  
   U.S. Patent **US7502723B1**, 2009.

3. **D. C. Jiles and D. L. Atherton**,  
   “Theory of Ferromagnetic Hysteresis,”  
   *Journal of Magnetism and Magnetic Materials*,  
   Vol. 61, 1986, pp. 48–60.

4. **LTspice**, Analog Devices.

---
