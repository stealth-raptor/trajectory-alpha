# UR5 PID and FOPID Project Guide

## Purpose

This project compares two joint-position controllers for a six-joint UR5 robotic arm model:

1. Conventional PID.
2. Fractional-order PID (FOPID).

FBPA optimization, standalone debug scripts, standalone test scripts and duplicate output scripts have been removed. The single entry point is:

```matlab
run_paper_outputs
```

On the first run, the script performs the bounded PID/FOPID gain search and
saves the result to `tuned_pid_fopid_gains.mat`. Later runs load that file and
reuse the same gains without repeating optimization. To intentionally tune
again, set `cfg.forceRetune = true` in `run_paper_outputs.m`, run once, then
set it back to `false`.

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
UR5 plant model
        |
        v
Position and velocity response
        |
        v
Metrics, plots and saved result files
```

The current plant prefers MATLAB's coupled `rigidBodyTree` `universalUR5` model. If Robotics System Toolbox is unavailable or the model fails its numerical preflight check, the plant uses a documented Table-2-based decoupled approximation.

## Running the project

From the project folder, run:

```matlab
run_paper_outputs
```

The script uses:

- Simulation time: 0 to 5 seconds.
- Sampling time: 0.002 seconds.
- Step: 30 degrees at 1 second.
- Sine: 30-degree amplitude and 0.5 Hz frequency.
- Six joints controlled independently at the controller level.
- RK4 numerical integration.
- Oustaloup approximation for the fractional operators.

The paper does not publish all experiment settings or final controller gains. The script keeps these assumptions visible in its configuration section.

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

The output includes:

- Six step-position tracking figures.
- Six step-error figures.
- One six-joint step-torque figure.
- Six sine-position tracking figures.
- Six sine-error figures.
- One six-joint sine-torque figure.
- CSV tables for joint-level and summary metrics.
- `paper_pid_fopid_results.mat` containing all simulation data.

## File descriptions

### `run_paper_outputs.m`

The only experiment script. It defines the simulation settings and baseline gains, creates step and sine references, runs PID and FOPID, calculates metrics, creates figures and saves all outputs.

It does not load or generate FBPA data.

### `make_ur5_plant.m`

Creates the robot plant.

When Robotics System Toolbox is available, it loads:

```matlab
loadrobot('universalUR5','DataFormat','column', ...
    'Gravity',[0 0 -9.81])
```

The rigid-body branch calculates acceleration using:

$$
\ddot q = M(q)^{-1}\left[\tau - h(q,\dot q) - G(q) - f(\dot q)\right]
$$

where:

- `M(q)` is the coupled inertia matrix.
- `h(q,dq)` contains velocity-product, Coriolis and centrifugal terms.
- `G(q)` is the gravity torque.
- `f(dq)` is the added friction torque.

The code prints the active model description. The preferred output is:

```text
coupled rigidBodyTree universalUR5
```

If the Robotics System Toolbox path is unavailable, the fallback uses the paper's listed masses, lengths and diagonal inertia values. The fallback is explicitly labelled:

```text
decoupled Table-2 approximation (not coupled UR5 dynamics)
```

### `simulate_ur5_controller.m`

Runs the six-joint closed-loop simulation. At each time step it:

1. Reads joint position and velocity.
2. Calculates the position error.
3. Calculates PID or FOPID torque.
4. Applies torque saturation.
5. Integrates the robot dynamics.
6. Stores position, velocity and torque.

Its result structure contains:

```matlab
result.t
result.q
result.dq
result.tau
```

### `init_oustaloup_bank.m`

Initializes the discrete filter banks used to approximate the fractional integral and derivative in FOPID. The default approximation uses order 5, which produces 11 first-order sections per filter.

### `oustaloup_step_bank.m`

Advances the Oustaloup filter states by one sample and returns the fractional integral and derivative signals used by FOPID.

### `make_step_metrics.m`

Calculates, for every joint and controller:

- Overshoot percentage.
- Settling time.
- Peak time.
- Target value.

It also creates a summary table averaged over the six joints.

### `make_sine_metrics.m`

Calculates, for every joint and controller:

- Mean squared error.
- Root mean squared error.
- Mean absolute error.
- Maximum absolute error.
- Mean absolute torque.
- Maximum absolute torque.

### `paper.txt`

Text copy of the research paper used as the reference for the controller comparison, UR5 parameters, performance metrics and experiment motivation.

## Controller explanation

### PID

The conventional PID control law is:

$$
\tau(t) = K_p e(t) + K_i\int e(t)dt + K_d\frac{de(t)}{dt}
$$

`Kp`, `Ki` and `Kd` are defined separately for each of the six joints.

### FOPID

FOPID adds two fractional orders:

$$
\tau(t) = K_p e(t) + K_i D^{-\lambda}e(t) + K_d D^{\mu}e(t)
$$

The five parameters for each joint are:

```text
Kp, Ki, Kd, lambda, mu
```

`lambda` controls the fractional integral order and `mu` controls the fractional derivative order. In this cleaned project, FOPID uses explicit baseline values; FBPA is no longer part of the workflow.

## Relationship to the research paper

The paper describes a Lagrange-based rigid-body model and gives some physical parameters. It also compares PID, FOPID and FBPA-FOPID. This cleaned project intentionally reproduces only the PID-versus-FOPID part.

The paper does not provide enough information to recreate its exact Simulink model or final gains. Therefore:

- The preferred plant is MATLAB's built-in coupled `universalUR5` model.
- The fallback uses the paper's Table 2 values but is decoupled and approximate.
- The reference amplitudes, frequency, friction, torque limit, initial conditions and gains are explicit implementation assumptions.
- The generated results should be presented as a MATLAB reproduction-style comparison, not as exact numerical reproduction of the paper.

## How to explain the project to a supervisor

> The project evaluates joint-space trajectory tracking for a six-DOF UR5 model. A conventional PID controller is used as the baseline and is compared with a fractional-order PID controller. The FOPID controller introduces fractional integral and derivative orders, giving two additional tuning parameters per joint. The robot is simulated using MATLAB's coupled UR5 rigid-body model when available, with a Table-2 approximation as a fallback. Step and sinusoidal references are applied, and the controllers are compared using overshoot, settling time, peak time, tracking error and torque effort.

The most important limitation is that the exact Simulink model and final gains from the paper are not published. The project therefore provides a transparent and reproducible comparison under stated assumptions.
