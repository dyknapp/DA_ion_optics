% =========================================================================
% Einzel Lens Simulation using Differential Algebra (DA) and DoPri8(7) Integration
% MODIFIED TO INCLUDE SEGMENTED QUADRUPOLE BENDER DA MAPS
% =========================================================================
% This script simulates particle trajectories through an Einzel lens system
% including a 90-degree quadrupole bender section.
% It uses:
% 1. BEM for electrostatic fields in straight sections.
% 2. Pre-calculated segmented linear maps for the quadrupole bender, converted to DA maps.
% 3. Differential Algebra (DA) to represent particle state and potential/maps.
% 4. Dormand-Prince 8(7) (DoPri8(7)) to integrate the DA map through straight sections.
% 5. Sequential evaluation of DA maps for ray tracing through the combined system,
%    including visualization within the bender.
% 6. STL generation for visualization (optional).
%
% Dependencies:
% - DiffAlg class, Electrode class, BEM_monopole, axial_potential, ...
% - run_dopri87_da_integration, integrate_dopri87_da, dopri87_step_da, get_dopri87_coeffs
% - trace_ray_fan_multisegment_da_map (New function)
% - generate_lens_trajectory_stl, stlwrite
% - bender_segmented_data.mat (Generated by modified quadrupole_bender_tracer.m)

set(0,'DefaultFigureWindowStyle','docked')
clear; clc;

% --- Basic Simulation Parameters ---
simParams.res = 1;        % Resolution for electrode boundary discretization
simParams.spacing = 0;    % Spacing between electrodes
simParams.z0 = 0;         % Starting z position for the first electrode
simParams.R = 1;          % Inner radius of the tube electrodes

% --- Differential Algebra & Integration Parameters ---
daParams.daOrder = 5;       % DA order (e.g., 3 for up to 3rd order effects)
daParams.nVars = 2;         % Number of variables (r, p_r)
daParams.q_charge = 1.0;    % Particle charge (e.g., +1.0 for proton)
daParams.m_mass = 1.0;      % Particle mass
daParams.E0 = 0.6;          % Initial Kinetic Energy (reference for potential energy U=qV-E0)

dopriParams.z_start_overall = 28.75;   % Overall start z for simulation
dopriParams.z_end_overall = 477.7;    % Overall end z for simulation
% Note: Step size (dz_step) is now handled within run_dopri87_da_integration options

% --- Define Integrator Options ---
integratorOpts = struct();
integratorOpts.RelTol = 1e-8;
integratorOpts.AbsTol = 1e-11;
integratorOpts.Stats = 'on'; % Get stats
integratorOpts.MaxStep = 2.0; % Example MaxStep


% Electrode geometry definition
BENDER_PLACEHOLDER = 0.5 * 43 * pi + 7.95*2;
tt = 2.0; pt = 4.0;
ids = {15.5,20,17,16,6,32};
ls  = {33.0,46,12,12,6,26};
ts  = {33,46+4.5,12+4.5,12+4.5,6+4.5,26};

electrodeData = { ...
struct('thickness',11.5,'shape',@(res)ETH_endcap(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_profile(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_profile(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_profile(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_profile(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_narrow(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_profile(res)), ...
struct('thickness',23.5,'shape',@(res)ashfold_steer_tube(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_profile(res)), ...
struct('thickness',11.5,'shape',@(res)ashfold_empty_tube(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_endcap_right(res)), ...
{'space', BENDER_PLACEHOLDER}, ...
struct('thickness',11.5,'shape',@(res)ETH_endcap(res)), ...
struct('thickness',23.5,'shape',@(res)ashfold_steer_tube(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_endcap_right(res)), ...
struct('thickness',55.0,'shape',@(res) tube(25, 50, 32*res)), ...
struct('thickness',ts{1},'shape',@(res)tube_in_plate(ids{1},tt,ls{1},pt,37.5,res)), ...
struct('thickness',ts{2},'shape',@(res)tube_in_plate(ids{2},tt,ls{2},pt,37.5,res)), ...
struct('thickness',ts{3},'shape',@(res)tube_in_plate(ids{3},tt,ls{3},pt,37.5,res)), ...
struct('thickness',ts{4},'shape',@(res)tube_in_plate(ids{4},tt,ls{4},pt,37.5,res)), ...
struct('thickness',ts{5},'shape',@(res)tube_in_plate(ids{5},tt,ls{5},pt,37.5,res)), ...
struct('thickness',ts{6},'shape',@(res)tube_in_plate(ids{6},tt,ls{6},pt,37.5,res)), ...
{'space',1}, ...
struct('thickness',8,'shape',@(res)tube_in_plate(3,0.5,6,2,37.5,res)) ...
};
% OPTIMIZATION PARAMS
el1 =  0.318;
se1 =  0.198;
se2 =  0.051;
dl1 = -1.428;
dl2 = -1.989;
dl3 = -1.654;
dl4 = -0.378;
ecv = -1.410;
bemParams.Vs = [1.0; 0.75; 0.50; 0.25; ...    % 1-4: uniform field tube
                0.0; el1; 0.0; ...            % 5-7: Einzel 1
                se1; 0.0; 0.0; 0.0; ...       % 8-11: steering electrode 1
                % No voltage for the bender 'space' (item 10)
                0.0; se2; 0.0;...             % 10-12: Einzel 2
                0.0; ...
                0.0; dl1; dl2; dl3; dl4; ecv; % Decelerator assembly
                -0.5; % endcap electrode
               ];

%% --- Load Bender Segmented Data ---
fprintf('Loading pre-calculated segmented bender data...\n');
bender_segmented_data_file = 'MAT007_clean_DA/bender_segmented_data.mat';
try
    loaded_data = load(bender_segmented_data_file, 'transfer_maps', 'ref_states', 'z_coords');
    transfer_maps_bender = loaded_data.transfer_maps;
    ref_states_bender = loaded_data.ref_states;
    z_coords_bender = loaded_data.z_coords;
    fprintf('  Loaded transfer_maps, ref_states, z_coords.\n');
    n_bender_segments = numel(transfer_maps_bender);
    n_bender_points = size(ref_states_bender, 2);
    n_bender_zcoords = numel(z_coords_bender);
    if n_bender_segments ~= (n_bender_points - 1) || n_bender_points ~= n_bender_zcoords
        error('Loaded segmented bender data has inconsistent dimensions.');
    end
    if ~iscell(transfer_maps_bender) || ~ismatrix(ref_states_bender) || ~isvector(z_coords_bender)
         error('Loaded segmented bender data has unexpected types.');
    end
    bender_length_loaded = z_coords_bender(end) - z_coords_bender(1);
    fprintf('  Bender length from loaded z_coords: %.6f\n', bender_length_loaded);
catch ME_load
    fprintf('Error loading %s: %s\n', bender_segmented_data_file, ME_load.message);
    error('Could not load necessary segmented bender data.');
end
fprintf('Segmented bender data loaded successfully.\n\n');

% --- Update Electrode Data with Loaded Bender Length & Create Stacks ---
fprintf('Updating electrodeData with loaded bender length...\n');
bender_idx = -1;
for idx = 1:length(electrodeData)
    if iscell(electrodeData{idx}) && strcmpi(electrodeData{idx}{1}, 'space')
        bender_idx = idx;
        % electrodeData{idx}{2} = bender_length_loaded; % Update length
        fprintf('  Updated electrodeData item %d (bender space) with loaded length %.6f\n', idx, bender_length_loaded);
        break;
    end
end
if bender_idx == -1
     error("Could not find the {'space', ...} element for the bender in electrodeData.");
end

% Define options for createElectrodeStack
opts_stack = struct('z0', simParams.z0, 'spacing', simParams.spacing, 'res', simParams.res);

% Build the FULL electrode stack (including space) for visualization and boundaries
fprintf('\nBuilding FULL electrode stack (including bender space)');
[electrodes_full_stack, z_boundaries_elements] = createElectrodeStack(electrodeData, opts_stack);
fprintf('Full electrode stack created (%d objects, %d boundaries).\n', numel(electrodes_full_stack), numel(z_boundaries_elements));

% Determine Bender Z-Range using the returned boundaries
fprintf('\nDetermining z-range for bender from calculated boundaries');
if bender_idx < 1 || bender_idx+1 > numel(z_boundaries_elements)
    error('Bender index is out of bounds for z_boundaries_elements.');
end
z_bend_start = 141.7;
z_bend_end = pi*86/4 + 7.95 + z_bend_start;
fprintf('  Bender insertion range: z = [%.3f, %.3f]\n\n', z_bend_start, z_bend_end);
z_coords_bender = z_coords_bender + z_bend_start;

% Check Voltage assignments against BEM electrodes
n_electrodes_structs_bem = numel(electrodes_full_stack);
if n_electrodes_structs_bem ~= length(bemParams.Vs)
     error("Voltage settings (%d) don't match the number of electrode structs for BEM (%d)!", length(bemParams.Vs), n_electrodes_structs_bem);
else
     fprintf('Voltage assignment count (%d) matches number of BEM electrodes (%d).\n\n', length(bemParams.Vs), n_electrodes_structs_bem);
end

% The loop calculating z_boundaries_elements is removed as it's now done by createElectrodeStack
% --- BEM Analysis (on stack excluding bender) ---
fprintf('Running BEM analysis...\n');
[qs, bemTable] = BEM_monopole(electrodes_full_stack, 1.0, false);
fprintf('BEM analysis complete.\n\n');

%% --- Voltage Configuration and Potential Calculation ---
% bemParams.Vs = bemParams.Vs; % Already defined
bemParams.qVs = qs' * bemParams.Vs; % Calculate charge distribution for this voltage setting

% Calculate and plot analytical axial potential
fprintf('Calculating axial potential...\n');
z_plot_min = min(bemTable.z_center)-simParams.spacing*2;
z_plot_max = max(bemTable.z_center)+simParams.spacing*2;
potParams.zs_plot = linspace(z_plot_min, z_plot_max, 1024);
potParams.U_analytical = axial_potential(potParams.zs_plot, bemParams.qVs, bemTable.r_center, bemTable.z_center);

figure(1); clf;
plot(potParams.zs_plot, potParams.U_analytical, '-k', 'LineWidth', 1.5);
title('Axial Potential (Excluding Bender Region)'); xlabel('Axial Position z'); ylabel('Potential U(z)');
grid on; hold on;
xline([z_bend_start, z_bend_end], '--r', {'Bender Start', 'Bender End'});
hold off; drawnow;
fprintf('Potential calculation and plotting complete.\n\n');

% ==========================================================
% --- DA Integration (Segmented) ---
% ==========================================================

sim_daOrder = daParams.daOrder; % Use the order defined earlier
sim_nVars = daParams.nVars;     % Use nVars defined earlier (should be 2)

% --- Run DA Integration 1 (Start to Bender) ---
fprintf('Starting DA DoPri8(7) Integration 1 (z=%.2f to z=%.2f)...\n', dopriParams.z_start_overall, z_bend_start);
dopriParams1 = struct('z_start', dopriParams.z_start_overall, 'z_end', z_bend_start);
[z_steps1, S_steps1, stats1] = ...
    run_dopri87_da_integration(dopriParams1, daParams, bemParams, bemTable, integratorOpts);
integrationResults1 = struct('z_steps', z_steps1, 'S_steps', {S_steps1}, 'stats', stats1);
fprintf('DA Integration 1 finished.\n');
fprintf('  Integration Stats: %d steps, %d failed, %d function evals.\n', ...
        stats1.nSteps, stats1.nFailed, stats1.nFunEvals);
if abs(z_steps1(end) - z_bend_start) > 1e-6
     warning('Integration 1 did not reach z_bend_start accurately (ended at %.6f).', z_steps1(end));
end

% --- Construct Bender DA Map Sequence ---
fprintf('\nConstructing DA Map Sequence for Bender Section...\n');
r_da = DiffAlg.var(1, sim_daOrder, sim_nVars);
pr_da = DiffAlg.var(2, sim_daOrder, sim_nVars);
one_da = DiffAlg.one(sim_nVars, sim_daOrder);

n_bender_segments = numel(transfer_maps_bender);
S_maps_bender_segmented = cell(1, n_bender_segments); % Cell array to hold DA maps for each segment

for i = 1:n_bender_segments
    M_i = transfer_maps_bender{i};  % 2x2 linear map for segment i
    X_ref_i = ref_states_bender(:, i); % Reference state [r; p] at START of segment i
    X_ref_f = ref_states_bender(:, i+1);% Reference state [r; p] at END of segment i

    % Calculate constant offset O_i = X_ref_f - M_i * X_ref_i
    Offset_i = X_ref_f - M_i * X_ref_i;
    Or_i = Offset_i(1);
    Op_i = Offset_i(2);

    % Construct DA map for segment i: X_f = M_i * X_i + O_i
    r_map_i = r_da * M_i(1,1) + pr_da * M_i(1,2) + one_da * Or_i;
    pr_map_i = pr_da * M_i(2,2) + r_da * M_i(2,1) + one_da * Op_i; % Reordered pr_da term

    % Store the DA map pair for this segment
    S_maps_bender_segmented{i} = {r_map_i, pr_map_i};

    % Optional: Check map at reference point (should recover X_ref_f)
    % r_check = r_map_i.evaluate(X_ref_i');
    % p_check = pr_map_i.evaluate(X_ref_i');
    % if abs(r_check - X_ref_f(1)) > 1e-9 || abs(p_check - X_ref_f(2)) > 1e-9
    %     warning('Bender DA map check failed for segment %d.', i);
    % end
end
fprintf('Bender DA map sequence constructed (%d segments, Order %d).\n', n_bender_segments, sim_daOrder);

% --- Run DA Integration 2 (Bender End to Overall End) ---
fprintf('\nStarting DA DoPri8(7) Integration 2 (z=%.2f to z=%.2f)...\n', z_bend_end, dopriParams.z_end_overall);
dopriParams2 = struct('z_start', z_bend_end, 'z_end', dopriParams.z_end_overall);
[z_steps2, S_steps2, stats2] = ...
    run_dopri87_da_integration(dopriParams2, daParams, bemParams, bemTable, integratorOpts);
integrationResults2 = struct('z_steps', z_steps2, 'S_steps', {S_steps2}, 'stats', stats2);
fprintf('DA Integration 2 finished.\n');
fprintf('  Integration Stats: %d steps, %d failed, %d function evals.\n', ...
        stats2.nSteps, stats2.nFailed, stats2.nFunEvals);
if abs(z_steps2(end) - dopriParams.z_end_overall) > 1e-6
     warning('Integration 2 did not reach z_end_overall accurately (ended at %.6f).', z_steps2(end));
end

%%
% ==========================================================
% RAY TRACING (Using Multi-Segment Maps)
% ==========================================================
fprintf('\n--- Multi-Segment Ray Tracing ---\n');

% --- Flexible Ray Tracing Parameters ---
traceParams.num_r_rays = 2;         % Number of rays in position (e.g., 5)
traceParams.num_s_rays = 3;         % Number of rays in slope (e.g., 5)
traceParams.max_r0_dev = 0.5;       % Max deviation in r from center_r0
traceParams.max_s0_dev = 0.01;      % Max deviation in p_r from center_s0
traceParams.center_r0 = 0.5;        % Central starting r position
traceParams.center_s0 = 0.0;        % Central starting p_r slope

% --- Specify DA Evaluation Order ---
trace_eval_order = sim_daOrder; % Use full order by default, can be reduced <= sim_daOrder

% --- Perform Ray Tracing using the new multi-segment function ---
fprintf('Performing ray tracing using multi-segment DA map evaluation (Order <= %d)...\n', trace_eval_order);


% --- Prepare Visualization Parameters ---
% Pass the full stack geometry (including space) for visualization
visParams.electrodes_full_stack = electrodes_full_stack; % Use the stack created earlier
visParams.simParams = simParams;   % Pass sim parameters (like R)
% --- Call the new multi-segment tracing function ---
figure(2); clf; subplot(5,1,1:3);
[traceResults.r_positions, traceResults.z_positions, traceResults.initial_states] = ...
    trace_ray_fan_multisegment_da_map(traceParams, ...
                                      integrationResults1, ...
                                      S_maps_bender_segmented, ... % Pass the sequence of bender maps
                                      z_coords_bender, ...        % Pass the corresponding z coordinates
                                      integrationResults2, ...
                                      visParams, ...              % Pass the updated visParams struct
                                      trace_eval_order);

fprintf('Ray tracing and plotting complete (Figure 2).\n\n');
drawnow

% Relative potential plot.
V0 = interp1(potParams.zs_plot, potParams.U_analytical, dopriParams.z_start_overall);
% Since V0 is very large compared to E0, we can say that we are almost
% independent of V0 for the overall trajectory.
xlabel('')
set(gca, 'XTickLabel', [])
ylim([-20 35])

hold on
lowerlims = xlim;
subplot(5, 1, 4:5)
plot(potParams.zs_plot, potParams.U_analytical / V0, '-k', 'LineWidth', 1.5)
xlim(lowerlims)
ylabel('Axial Potential V/V0')
yline(0, '--k')
xline(dopriParams.z_start_overall, '--k', 'V0')
hold off
grid
xlabel('Position Along Beam Axis [mm]')
drawnow


%%
figure(3); clf;
subplot(1, 2, 1)
% --- Flexible Ray Tracing Parameters ---
traceParams.num_r_rays = 2;         % Number of rays in position (e.g., 5)
traceParams.num_s_rays = 3;         % Number of rays in slope (e.g., 5)
traceParams.max_r0_dev = 0.5;       % Max deviation in r from center_r0
traceParams.max_s0_dev = 0.005;      % Max deviation in p_r from center_s0
traceParams.center_r0 = 0.5;        % Central starting r position
traceParams.center_s0 = 0.0;        % Central starting p_r slope
trace_eval_order = 1;
[traceResults.r_positions, traceResults.z_positions, traceResults.initial_states] = ...
    trace_ray_fan_multisegment_da_map(traceParams, ...
                                      integrationResults1, ...
                                      S_maps_bender_segmented, ... % Pass the sequence of bender maps
                                      z_coords_bender, ...        % Pass the corresponding z coordinates
                                      integrationResults2, ...
                                      visParams, ...              % Pass the updated visParams struct
                                      trace_eval_order);
xlim([430 477.7]); ylim([-4,4])
title('Paraxial Result (Aberrations Omitted)')
text(467.5, -3.8, 'Ion Trap Endcap', 'Rotation', 90)
text(443,  3.8, 'Objective Lens', 'Rotation', 270)
text(432,  0.75, 'r_0 = 0 mm', 'Color', [0.8, 0.2, 1.0])
text(432, -2.50, 'r_0 = 1 mm', 'Color', [0.2, 0.8, 1.0])

subplot(1, 2, 2)
% --- Flexible Ray Tracing Parameters ---
traceParams.num_r_rays = 2;         % Number of rays in position (e.g., 5)
traceParams.num_s_rays = 3;         % Number of rays in slope (e.g., 5)
traceParams.max_r0_dev = 0.5;       % Max deviation in r from center_r0
traceParams.max_s0_dev = 0.005;      % Max deviation in p_r from center_s0
traceParams.center_r0 = 0.5;        % Central starting r position
traceParams.center_s0 = 0.0;        % Central starting p_r slope
trace_eval_order = 5;
[traceResults.r_positions, traceResults.z_positions, traceResults.initial_states] = ...
    trace_ray_fan_multisegment_da_map(traceParams, ...
                                      integrationResults1, ...
                                      S_maps_bender_segmented, ... % Pass the sequence of bender maps
                                      z_coords_bender, ...        % Pass the corresponding z coordinates
                                      integrationResults2, ...
                                      visParams, ...              % Pass the updated visParams struct
                                      trace_eval_order);
xlim([430 477.7]); ylim([-4,4])
title('Full Result (Aberrations Included)')
print2pdf('figures/beam_plot_paraxial_detail.pdf', [15 10])
text(467.5, -3.8, 'Ion Trap Endcap', 'Rotation', 90)
text(443,  3.8, 'Objective Lens', 'Rotation', 270)
text(432,  0.75, 'r_0 = 0 mm', 'Color', [0.8, 0.2, 1.0])
text(432, -2.50, 'r_0 = 1 mm', 'Color', [0.2, 0.8, 1.0])

print2pdf('figures/beam_plot_paraxial_detail.pdf', [15, 10])

%%
% ==========================================================
% STL TRAJECTORY TUBE GENERATION
% ==========================================================
fprintf('\n--- Generating STL Visualization for Trajectories ---\n');

% --- Define Parameters for STL Generation ---
stlParams_traj.tube_radius = 0.1; % Adjust thickness of tubes (relative to R=1)
stlParams_traj.n_tube_segments = 12; % Smoothness of tubes (e.g., 8-16)
stlParams_traj.axis_z_min = -100; % Extend axis before simulation start
stlParams_traj.axis_z_max = 600;  % Extend axis after simulation end
stlParams_traj.output_folder = 'trajectory_tube_stl'; % Subfolder for output
stlParams_traj.output_filename = fullfile(stlParams_traj.output_folder, 'trajectory_fan_and_axis.stl');

% Define R0 used for the bender simulation (ensure this matches quadrupole_bender_tracer.m)
% This value needs to be explicitly set or loaded if not already available.
R0_bender = 86; % Example value, **ADJUST THIS TO MATCH YOUR BENDER SIM**

% --- Check if required data exists ---
if exist('traceResults', 'var') && isfield(traceResults, 'r_positions') && isfield(traceResults, 'z_positions') && ...
   exist('z_bend_start', 'var') && exist('z_bend_end', 'var')

    % --- Call the STL Generation Function ---
    try
        generate_trajectory_tubes_stl(traceResults.z_positions, ... % N x 1 z-coords
                                      traceResults.r_positions, ... % N x M r-coords
                                      z_bend_start, ...
                                      z_bend_end, ...
                                      R0_bender, ...
                                      stlParams_traj);
        fprintf('STL generation call completed.\n');
    catch ME_stl
        fprintf('!!! Error calling generate_trajectory_tubes_stl: %s\n', ME_stl.message);
        fprintf('    Location: %s, Line: %d\n', ME_stl.stack(1).name, ME_stl.stack(1).line);
    end
else
    warning('Required trajectory data (traceResults) or bender boundaries (z_bend_start/end) not found. Skipping trajectory STL generation.');
end

disp('--- Main Script Finished ---');
