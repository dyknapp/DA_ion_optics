% --- Setup ---
set(0,'DefaultFigureWindowStyle','docked') % Dock figures for easier management
clear; clc; clf;

fprintf('--- Minimal BEM Monopole & Quadrupole Script ---\n');

% --- Parameters ---
res = 4; % Resolution for electrode boundary discretization (lower for faster testing)

% --- Electrode Geometry Definition ---
fprintf('Defining electrode geometry...\n');
electrodeData = { ...
    struct('thickness', 14.2, 'shape', @(res) tube_in_plate(1.5, 6.5, 14, 9, 60, res)), ...
    {'space', -7.9}, ...
    struct('thickness', 50.0, 'shape', @(res) tube(11.8, 50, res*32)), ...
    {'space',  3.1}, ...
    struct('thickness', 14.2, 'shape', @(res) tube_in_plate(1.5, 6.5, 14, 9, 60, res)), ...
};

% --- Create Electrode Stack ---
% Use createElectrodeStack to handle geometry creation and positioning
fprintf('Building electrode stack...\n');
opts_stack = struct('z0', -31, 'spacing', 0, 'res', res);
try
    [electrodes, z_boundaries_elements] = createElectrodeStack(electrodeData, opts_stack);
    nElectrodes = numel(electrodes);
    fprintf('Electrode stack created with %d electrodes.\n', nElectrodes);
    fprintf('Element boundaries (z): %s\n', mat2str(z_boundaries_elements, 4));
catch ME
    fprintf('Error creating electrode stack: %s\n', ME.message);
    error('Failed to create electrode stack. Ensure Electrode.m and createElectrodeStack.m are available.');
end

% --- Define Voltages ---
Vs_monopole_open  = [-5.0;  5.0; 10.0];
Vs_monopole_close = [ 10.0; 5.0; 10.0];
Vs_quadpole = [  0.0; 450.0;  0.0]; % Column vector for quadrupole potentials

Vs_monopole_open = Vs_monopole_open(:); % Ensure column vector
Vs_quadpole = Vs_quadpole(:); % Ensure column vector

% --- Run BEM Monopole ---
fprintf('\nRunning Monopole BEM analysis...\n');
full_integration_cutoff_mono = 1.0;
verbose_bem = true; % Show BEM progress and results
try
    [qs_mono, bemTable_mono] = BEM_monopole(electrodes, full_integration_cutoff_mono, verbose_bem);
    fprintf('Monopole BEM analysis complete.\n');
catch ME_mono
    fprintf('Error during Monopole BEM: %s\n', ME_mono.message);
    error('Monopole BEM failed. Ensure BEM_monopole.m and its dependencies are available.');
end

% --- Run BEM Quadrupole ---
fprintf('\nRunning Quadrupole BEM analysis...\n');
full_integration_cutoff_quad = 1.0; % Can be same or different cutoff
try
    % We only need qs_quad, bemTable_quad will be similar to bemTable_mono
    % but with different charge columns. We use bemTable_mono for geometry.
    [qs_quad, ~] = BEM_quadrupole(electrodes, full_integration_cutoff_quad, verbose_bem);
    fprintf('Quadrupole BEM analysis complete.\n');
catch ME_quad
    fprintf('Error during Quadrupole BEM: %s\n', ME_quad.message);
    error('Quadrupole BEM failed. Ensure BEM_quadrupole.m and its dependencies are available.');
end

% --- Calculate Monopole Charge Distribution for Axial Potential ---
fprintf('\nCalculating monopole charge distribution for specified voltages...\n');
% qVs = sum over electrodes (qs_mono(elec,:) * Vs_monopole_open(elec))
% Or using matrix multiplication:
qVs_mono_open = qs_mono' * Vs_monopole_open; % Result is [numSegments x 1]
qVs_mono_close = qs_mono' * Vs_monopole_close; % Result is [numSegments x 1]

% --- Calculate Axial Potential (using Monopole results) ---
fprintf('Calculating axial potential (from monopole BEM)...\n');
% Define z-range for plotting potential
zs_plot = linspace(-35, 35, 1024); % Points for plotting

% Calculate potential using the axial_potential function
try
    U_axial_mono_open = axial_potential(zs_plot, qVs_mono_open, bemTable_mono.r_center, bemTable_mono.z_center);
    U_axial_mono_close = axial_potential(zs_plot, qVs_mono_close, bemTable_mono.r_center, bemTable_mono.z_center);
catch ME_axial
     fprintf('Error during axial potential calculation: %s\n', ME_axial.message);
     error('Axial potential calculation failed. Ensure axial_potential.m is available.');
end
fprintf('Axial potential calculated.\n');

%% --- Plot Axial Potential ---
figure(1); clf; % Create/clear figure 1
plot(zs_plot, U_axial_mono_open, '-b', 'LineWidth', 1.5);
hold on;

plot(zs_plot, U_axial_mono_close, '-.r', 'LineWidth', 1.5);
% plot(bemTable_mono.z_center, bemTable_mono.r_center, '.k')

% Add labels and title
xlabel('Axial Position [mm]');
ylabel('Axial Potential qU(z) [eV]');
title('Axial Potential of Ion Trap (Floated at +5V)');
grid on
xline(-23.5, '--k', 'LineWidth', 1, 'HandleVisibility', 'off')
xline( 23.5, '--k', 'LineWidth', 1)
yline(0, '-k', 'LineWidth', 1, 'HandleVisibility', 'off')
yline(5, '--b', 'LineWidth', 1)
hold off
legend('Trap Open', 'Trap Closed', 'Endcap Surface', 'Trap Offset', 'Location', 'best', 'FontSize', 16);

% fontsize('scale', 2)
