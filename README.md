# Flexible Link Simulation

MATLAB simulation code for **Vibration Suppression and Tracking Performance of a Rotary Flexible Link Manipulator under a Common Motion Deadline**.

## Requirements

MATLAB with Control System Toolbox. Original results used MATLAB R2025b Update 3.

## Run and save results

Open this folder in MATLAB and enter:

```matlab
deadline_simulation
```

This computes five configurations across 27 design cases, 48 evaluation cases, and six sensitivity scenarios. Parameters, full precision coefficients, and cases are included in the code. Results are saved in `results`:

- `ieee_model_study.mat`: time, reference, shaped commands, all 75 × 5 response histories, nominal responses, parameters, coefficients, metrics, and deadline errors.
- `ieee_stress_sensitivity.mat`: all 6 × 5 sensitivity response histories, scenario definitions, parameters, commands, metrics, and deadline errors.
- CSV files: nominal, evaluation, and sensitivity metrics and summaries.

In both MAT files, `responses{case,method}` contains columns `[theta, beta, Vref]`, with units `[rad, rad, V]`. Time and the original reference are in `meta.t` and `meta.reference`. Method order is in `meta.methods`. `commands(:,method)` stores the shaped reference in rad. The main study deadline errors are in rad; sensitivity deadline errors are in degrees. Rerunning overwrites matching output files.

## Plot saved results

```matlab
plot_results
```

This reads the saved results and exports Figs. 4–7 as PNG, PDF, and MATLAB FIG. It does not rerun the simulation.

## Optional redesign

```matlab
deadline_simulation('redesign')
```

Requires Robust Control Toolbox and Optimization Toolbox in addition. The default command uses the fixed manuscript coefficients. The code label `scenario LP` denotes the scenario optimized FIR shaper.

The revised files have not yet been executed in MATLAB after restructuring.
