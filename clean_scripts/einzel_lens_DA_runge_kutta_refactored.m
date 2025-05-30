% =========================================================================
% Einzel Lens Simulation using Differential Algebra (DA) and RK4 Integration
% =========================================================================
% This script simulates particle trajectories through an Einzel lens.
% It uses:
% 1. Boundary Element Method (BEM) to calculate the electric field.
% 2. Differential Algebra (DA) to represent particle state and potential.
% 3. 4th-Order Runge-Kutta (RK4) to integrate the DA map through the lens.
% 4. Ray tracing using the calculated transfer matrices.
% 5. STL generation for visualization.
%
% Dependencies:
% - DiffAlg class
% - Electrode class
% - BEM_monopole function
% - axial_potential function
% - axial_potential_coeffs function
% - createElectrodeStack function
% - run_rk4_da_integration function (New)
% - trace_ray_fan function (New)
% - generate_lens_trajectory_stl function (New)
% - stlwrite function (from File Exchange or built-in)

set(0,'DefaultFigureWindowStyle','docked')
clear; clc;

% --- Basic Simulation Parameters ---
simParams.res = 256;      % Resolution for electrode boundary discretization
simParams.spacing = 1;    % Spacing between electrodes
simParams.z0 = 0;         % Starting z position for the first electrode
simParams.R = 1;          % Inner radius of the tube electrodes

% --- Electrode Geometry Definition ---
% Define electrode shapes and thicknesses
electrodeData = {
    struct('thickness', 10, 'shape', @(res) tube(simParams.R, 10, res))
    struct('thickness', 10, 'shape', @(res) tube(simParams.R, 10, res))
    struct('thickness', 10, 'shape', @(res) tube(simParams.R, 10, res))
};

% Build the electrode stack
fprintf('Building electrode stack...\n');
electrodes = createElectrodeStack(electrodeData, 'z0', simParams.z0, 'spacing', simParams.spacing, 'res', simParams.res);
fprintf('Electrode stack created.\n\n');

% --- BEM Analysis ---
fprintf('Running BEM analysis...\n');
[qs, bemTable] = BEM_monopole(electrodes, 1.0, true); % Use verbose BEM output
fprintf('BEM analysis complete.\n\n');

% --- Voltage Configuration and Potential Calculation ---
bemParams.Vs = [0; -0.5; -1.0]; % Voltages applied to electrodes [V1; V2; V3]
bemParams.qVs = qs' * bemParams.Vs; % Calculate charge distribution for this voltage setting

% Calculate and plot analytical axial potential
fprintf('Calculating axial potential...\n');
potParams.zs_plot = linspace(min(bemTable.z_center)-simParams.spacing, max(bemTable.z_center)+simParams.spacing, 1024);
potParams.U_analytical = axial_potential(potParams.zs_plot, bemParams.qVs, bemTable.r_center, bemTable.z_center);

figure(1); clf;
plot(potParams.zs_plot, potParams.U_analytical, '-k', 'LineWidth', 1.5);
title('Axial Potential and Series Approximation'); xlabel('Axial Position z'); ylabel('Potential U(z)');
grid on; hold on;

% --- Potential Series Approximation (for context/verification) ---
potParams.Nmax_coeffs_plot = 10; % Order for plotting series approx
potParams.z0_expansion_plot = 12.1; % Center for series approx plot
potParams.expansion_range_plot = 1.0; % Range around center for plot

Coeffs_plot = axial_potential_coeffs(potParams.Nmax_coeffs_plot, potParams.z0_expansion_plot, bemParams.qVs, bemTable.r_center, bemTable.z_center);

% Display coefficients (optional)
fprintf('Axial Potential Coefficients C_n around z=%.2f:\n', potParams.z0_expansion_plot);
for n = 0:potParams.Nmax_coeffs_plot
    fprintf('  C_%2d = % .6e\n', n, Coeffs_plot(n+1));
end

local_idcs = (potParams.zs_plot >= potParams.z0_expansion_plot - potParams.expansion_range_plot) & ...
             (potParams.zs_plot <= potParams.z0_expansion_plot + potParams.expansion_range_plot);
z_local = potParams.zs_plot(local_idcs);
delta_z = z_local - potParams.z0_expansion_plot;
U_approx_plot = zeros(size(z_local));
for n = 0:potParams.Nmax_coeffs_plot
    U_approx_plot = U_approx_plot + Coeffs_plot(n+1) * delta_z.^n;
end
plot(z_local, U_approx_plot, '--r', 'LineWidth', 1.5);
xline(potParams.z0_expansion_plot, '-.r', 'Label', 'Expansion Center');
hold off;
legend('BEM Potential', 'Taylor Series Approx.', 'Expansion Center');
fprintf('Potential calculation and plotting complete.\n\n');

% --- Differential Algebra RK4 Integration Parameters ---
daParams.daOrder = 3;       % DA order (e.g., 3 for up to 3rd order effects)
daParams.nVars = 2;         % Number of variables (r, p_r)
daParams.q_charge = 1.0;    % Particle charge (e.g., +1.0 for proton)
daParams.m_mass = 1.0;      % Particle mass
daParams.E0 = 1.0;          % Initial Kinetic Energy (reference for potential energy U=qV-E0)

rk4Params.z_start = 3;      % Start z for integration
rk4Params.z_end   = 28;     % End z for integration
rk4Params.dz_step = 0.1;    % Step size for RK4 integration slices

% --- Run RK4 DA Integration ---
fprintf('Starting DA RK4 integration...\n');
[rk4Results.all_R_matrices, rk4Results.z_boundaries, rk4Results.R_total] = ...
    run_rk4_da_integration(rk4Params, daParams, bemParams, bemTable);
fprintf('DA RK4 integration finished.\n\n');

% --- Ray Tracing Parameters ---
traceParams.num_rays = 11;              % Number of rays to trace
traceParams.max_r0 = simParams.R * 0.75; % Max initial radius relative to tube radius
traceParams.initial_p0 = 0.01;             % Initial momentum (slope dr/dz)

% --- Perform Ray Tracing and Plotting ---
fprintf('Performing ray tracing...\n');
visParams.electrodes = electrodes; % Pass electrode geometry for plotting
visParams.simParams = simParams;   % Pass sim parameters (like R) if needed for plot limits

[traceResults.r_positions] = trace_ray_fan(traceParams, rk4Results, visParams);
fprintf('Ray tracing and plotting complete (Figure 2).\n\n');

drawnow

%% STL

% --- STL Generation Parameters ---
stlParams.R_outer = 5;                 % Outer radius for electrode block visualization
stlParams.num_angular_steps = 64;      % Rotational smoothness
stlParams.cutaway_angle_degrees = 0;  % Cutaway angle (0 for full revolution)
stlParams.output_folder = 'electrode_and_traj_meshes_stl'; % Output folder
stlParams.export_trajectory = true;    % Flag to include trajectory in STL
stlParams.export_lens = false;

% --- Generate STL Geometry ---
fprintf('Generating STL geometry...\n');
generate_lens_trajectory_stl(stlParams, visParams.electrodes, rk4Results.z_boundaries, traceResults.r_positions);
fprintf('STL generation complete.\n\n');

disp('--- Main Script Finished ---');
