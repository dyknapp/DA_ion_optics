set(0,'DefaultFigureWindowStyle','docked')
clear; clc;

% --- Basic Simulation Parameters ---
simParams.res = 2;        % Resolution for electrode boundary discretization
simParams.spacing = 0;    % Spacing between electrodes
simParams.z0 = 0;         % Starting z position for the first electrode
simParams.R = 1;          % Inner radius of the tube electrodes

% --- Differential Algebra & Integration Parameters ---
daParams.daOrder = 5;       % DA order (e.g., 3 for up to 3rd order effects)
daParams.nVars = 2;         % Number of variables (r, p_r)
daParams.q_charge = 1.0;    % Particle charge (e.g., +1.0 for proton)
daParams.m_mass = 1.0;      % Particle mass
daParams.E0 = 0.1;          % Initial Kinetic Energy (reference for potential energy U=qV-E0)

dopriParams.z_start_overall = 23.25;   % Overall start z for simulation
dopriParams.z_end_overall = 480.0;    % Overall end z for simulation
% Note: Step size (dz_step) is now handled within run_dopri87_da_integration options

% --- Define Integrator Options ---
integratorOpts = struct();
integratorOpts.RelTol = 1e-10;
integratorOpts.AbsTol = 1e-11;
integratorOpts.Stats = 'on'; % Get stats
integratorOpts.MaxStep = 1.0;

% Define Electrode Stack

electrodeData = { ...
    struct('thickness', 24, 'shape', @(res) ashfold_repeller_profile(res)), ...
    struct('thickness', 24, 'shape', @(res) ashfold_extractor_profile(res)), ...
    struct('thickness', 12, 'shape', @(res) ashfold_empty_tube(res)), ...
    struct('thickness', 12, 'shape', @(res) ashfold_plate_profile(res)), ...
    struct('thickness', 12, 'shape', @(res) ashfold_empty_tube(res)), ...
    struct('thickness',360, 'shape', @(res) tube(25, 395, 32*res)), ...
};

bemParams.Vs = [ ...
        2000.0; ...
        1650.0; ...
        0.0; ...
        -1000.0; ...
        0.0; ...
        0.0; ...
    ];

%% Set up the BEM calculation
% Define options for createElectrodeStack
opts_stack = struct('z0', simParams.z0, 'spacing', simParams.spacing, 'res', simParams.res);

% Build the FULL electrode stack (including space) for visualization and boundaries
fprintf('\nBuilding FULL electrode stack (including bender space)');
[electrodes, z_boundaries_elements] = createElectrodeStack(electrodeData, opts_stack);
fprintf('Full electrode stack created (%d objects, %d boundaries).\n', numel(electrodes), numel(z_boundaries_elements));

% Check Voltage assignments against BEM electrodes
n_electrodes_structs_bem = numel(electrodes);
if n_electrodes_structs_bem ~= length(bemParams.Vs)
     error("Voltage settings (%d) don't match the number of electrode structs for BEM (%d)!", length(bemParams.Vs), n_electrodes_structs_bem);
else
     fprintf('Voltage assignment count (%d) matches number of BEM electrodes (%d).\n\n', length(bemParams.Vs), n_electrodes_structs_bem);
end

%% Run BEM
fprintf('Running BEM analysis...\n');
[qs, bemTable] = BEM_monopole(electrodes, 1.0, false);
fprintf('BEM analysis complete.\n\n');

%% Prep the potential
% Calculate charge distribution for the specific voltage settings
bemParams.qVs = qs' * bemParams.Vs;

% Account for initial electrostatic potential
V0 = axial_potential_coeffs(0, simParams.z0, bemParams.qVs, ...
    bemTable.r_center, bemTable.z_center);
daParams.E0 = daParams.E0 + V0;
fprintf('\nInitial potential energy V0 = %5e\n\n', V0)

%% DA integration
fprintf('Starting DA DoPri8(7) Integration 1 (z=%.2f to z=%.2f)...\n', ...
    dopriParams.z_start_overall, dopriParams.z_end_overall);
dopriParams1 = struct('z_start', dopriParams.z_start_overall, 'z_end', dopriParams.z_end_overall);

% Integration
[z_steps, S_steps, stats] = ...
    run_dopri87_da_integration(dopriParams1, daParams, bemParams, bemTable, integratorOpts);
integrationResults = struct('z_steps', z_steps, 'S_steps', {S_steps}, 'stats', stats);

fprintf('DA Integration 1 finished.\n');
fprintf('  Integration Stats: %d steps, %d failed, %d function evals.\n', ...
        stats.nSteps, stats.nFailed, stats.nFunEvals);

%% Ray Trace using DA transfer maps
% --- Flexible Ray Tracing Parameters ---
traceParams.num_r_rays = 3;         % Number of rays in position (e.g., 5)
traceParams.num_s_rays = 3;         % Number of rays in slope (e.g., 5)
traceParams.max_r0_dev = 5;       % Max deviation in r from center_r0
traceParams.max_s0_dev = 0.1;      % Max deviation in p_r from center_s0
traceParams.center_r0 = 0.0;        % Central starting r position
traceParams.center_s0 = 0.0;        % Central starting p_r slope

% --- Specify DA Evaluation Order ---
trace_eval_order = 1; % Plot aberration-free

% --- Perform Ray Tracing using the new multi-segment function ---
fprintf('Performing ray tracing using multi-segment DA map evaluation (Order <= %d)...\n', trace_eval_order);


% --- Prepare Visualization Parameters ---
% Pass the full stack geometry (including space) for visualization
visParams.electrodes = electrodes;
visParams.simParams = simParams;   % Pass sim parameters (like R)
% --- Call the ray tracing function ---
figure(2); clf;
trace_ray_fan_da_map(traceParams, integrationResults, visParams, trace_eval_order);
fprintf('Ray tracing and plotting complete (Figure 2).\n\n');

%% Find focal length

