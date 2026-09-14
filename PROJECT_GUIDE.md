# UR5 PID and FOPID Project Guide

## Purpose

This project compares two joint-position controllers for a six-joint UR5 robotic arm model:

1. Conventional PID.
2. Fractional-order PID (FOPID).

The single entry point is:

```matlab
run_paper_outputs
```

The PID and FOPID gains are hardcoded directly in `run_paper_outputs.m`. The
project does not run an optimizer; changing the gains is an intentional manual
experiment.

This generates the paper-style step and sine tracking results for PID and FOPID.

## Main workflow

```text
Reference trajectory
        |
        v
Tracking error: e = q_ref - q
        |
        v
PID or FOPID controller
        |
        v
Joint torque tau
        |
        v
UR5 plant model  (paper-faithful Lagrange, fully coupled)
        |
        v
Position and velocity response
        |
        v
Metrics, plots and saved result files
```

## Plant model (FINAL, single architecture)

The plant is a **paper-faithful, fully-coupled Lagrange model**.

Dynamic equation (exactly as required by the reference paper):

```
M(q) * qdd + C(q,qd)*qd + G(q) + tau_f(qd) = tau
```

with

```
tau_f,i = fc_i * sign(qd_i) + b_i * qd_i
```

- All masses, lengths and inertia matrices are taken **exactly** from Table 2 of the paper.
- The kinematic chain (DH parameters) and COM locations are **not** fully specified by the paper; a clearly documented standard UR5-like convention is used. See the header comments of `paper_lagrange_ur5.m` for the complete list of geometric assumptions.
- There is **no** `loadrobot('universalUR5')`, no `rigidBodyTree`, and no decoupled single-joint approximation.

The plant is created by:

```matlab
plant = make_ur5_plant(cfg);
```

and exposes the acceleration map expected by the simulator:

```matlab
qdd = plant.accel(q, dq, tau);
```

## Running the project

From the project folder, run:

```matlab
run_paper_outputs
```

Validation of the plant can be performed with:

```matlab
validate_paper_plant
```

The script uses:

- Simulation time: 0 to 5 seconds.
- Sampling time: 0.002 seconds.
- Step: 30 degrees at 1 second.
- Sine: 30-degree amplitude and 0.5 Hz frequency.
- Six joints controlled independently at the controller level.
- RK4 numerical integration.
- Oustaloup approximation for the fractional operators.

## Output files

Results are written to:

```text
outputs/paper/
├── step/
│   ├── figures/
│   └── metrics/
├── sine/
│   ├── figures/
│   └── metrics/
└── results/
```

## File descriptions

### `run_paper_outputs.m`

The only experiment script. It defines the simulation settings and baseline gains, creates step and sine references, runs PID and FOPID, calculates metrics, creates figures and saves all outputs.

### `make_ur5_plant.m`

Creates the single definitive plant. Calls the Lagrange dynamics implementation and returns a struct with the `accel` handle.

### `paper_lagrange_ur5.m`

Core dynamics. Computes the configuration-dependent inertia matrix \(M(q)\), the Coriolis/centrifugal vector \(C(q,\dot q)\dot q\), the gravity vector \(G(q)\) and the friction vector from the Lagrange formulation (homogeneous transforms + geometric Jacobians of the centres of mass). All paper-supplied numerical data and all geometric assumptions are isolated inside this file.

### `simulate_ur5_controller.m`

Unchanged. Runs the six-joint closed-loop simulation using the plant's `accel` map.

### `validate_paper_plant.m`

Numerical verification that \(M\) is 6×6, symmetric and positive-definite, that the model is coupled, that \(G\) and friction have the correct structure, and that no toolbox dependency remains.

### Other files

`init_oustaloup_bank.m`, `oustaloup_step_bank.m`, `make_step_metrics.m`, `make_sine_metrics.m` are unchanged.

## Controller explanation

### PID

```
tau(t) = Kp * e(t) + Ki * integral(e) + Kd * de/dt
```

### FOPID

```
tau(t) = Kp * e(t) + Ki * D^(-lambda) e(t) + Kd * D^(mu) e(t)
```

Both controllers use model-based gravity compensation in the paper-style
experiments. The feedforward term `G(q)` is added to the PID/FOPID torque
before torque saturation, so the static gravity load does not need to be
rejected by the feedback gains alone. Set `cfg.useGravityCompensation` to
`false` to reproduce the uncompensated behavior.

## Relationship to the research paper

The paper describes a Lagrange-based rigid-body model and gives the physical parameters of Table 2. It does **not** publish a complete DH table, COM locations or the numerical friction coefficients. Consequently:

- Masses, lengths and the six diagonal inertia matrices are used exactly as published.
- Geometric quantities required by the Lagrange derivation but absent from the paper are supplied by a single, clearly documented set of assumptions (see `paper_lagrange_ur5.m`).
- The resulting plant is the definitive model for all subsequent PID / FOPID experiments in this repository.

The generated results should be presented as a MATLAB reproduction under these explicit modelling choices, not as a bit-exact numerical replica of the authors' unpublished Simulink model.
