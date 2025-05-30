% test_uniform_field_map.m
%
% Tests the RK4 DA integration by calculating the overall transfer map
% for a uniform electric field E_z and comparing its predictions for a
% fan of initial conditions against the analytical solution.

set(0,'DefaultFigureWindowStyle','docked')
clear; clc; close all;

% --- Simulation Parameters ---
% Physical Constants
params.q_charge = 1.0;      % Particle charge
params.m_mass   = 1.0;      % Particle mass
params.Ez       = 1.0;     % Uniform electric field strength along z (negative = accelerating for q>0)

% Initial Conditions for the REFERENCE trajectory
ref_init.r0  = 0.0;
ref_init.z0  = 0.0;
ref_init.vr0 = 0.1;  % Initial radial velocity
ref_init.vz0 = 1.0;  % Initial axial velocity

% Calculate derived initial conditions for the REFERENCE trajectory
ref_init.E0   = 0.5 * params.m_mass * (ref_init.vr0^2 + ref_init.vz0^2) - params.q_charge * params.Ez * ref_init.z0;
ref_init.pr0  = ref_init.vr0 / ref_init.vz0; % Initial pr = dr/dz

fprintf('--- Uniform Field Map Test Parameters ---\n');
fprintf('q=%.2f, m=%.2f, Ez=%.2f\n', params.q_charge, params.m_mass, params.Ez);
fprintf('Reference Initial State: r0=%.3f, z0=%.3f, vr0=%.3f, vz0=%.3f\n', ref_init.r0, ref_init.z0, ref_init.vr0, ref_init.vz0);
fprintf('Reference E0 = %.4f, Reference pr0 = %.4f\n', ref_init.E0, ref_init.pr0);

% RK4 Parameters
rk4Params.z_start = ref_init.z0;
rk4Params.z_end   = 50.0;   % Choose an appropriate end point
rk4Params.dz_step = 0.05; % Integration step size (smaller step -> more accurate map)

% DA Parameters
daParams.daOrder = 7;      % Order for DA calculations (higher order -> more accurate map)
daParams.nVars   = 2;      % State variables: r, pr

% --- Check for immediate reflection for reference trajectory ---
U_start_ref = -params.q_charge * params.Ez * rk4Params.z_start - ref_init.E0;
if U_start_ref >= 0
    error('Reference initial condition leads to U >= 0 (%.3e). Particle cannot start, check parameters.', U_start_ref);
end

% --- Slice Definition ---
total_length = rk4Params.z_end - rk4Params.z_start;
if total_length <= 0, error('Total z range is non-positive.'); end
num_full_slices = floor(total_length / rk4Params.dz_step);
dz_remainder = total_length - num_full_slices * rk4Params.dz_step;
remainder_tolerance = 1e-9 * rk4Params.dz_step;
has_remainder = dz_remainder > remainder_tolerance;
num_total_slices = num_full_slices + has_remainder;

fprintf('RK4 Setup: z=[%.2f, %.2f], dz=%.4f, Slices=%d, DA Order=%d\n', rk4Params.z_start, rk4Params.z_end, rk4Params.dz_step, num_total_slices, daParams.daOrder);
if num_total_slices == 0, error('No slices defined.'); end

% --- Initialization for Map Calculation ---
r_da  = DiffAlg.var(1, daParams.daOrder, daParams.nVars); % r (DA var 1)
pr_da = DiffAlg.var(2, daParams.daOrder, daParams.nVars); % p_r = dr/dz (DA var 2)
S_identity = {r_da, pr_da}; % Initial state = Identity map

% --- Calculate Overall Transfer Map using DoPri8(7) ---
fprintf('Calculating overall DA transfer map using DoPri8(7)...\n');
tic;

% Initial state = Identity map
S_identity = {r_da, pr_da};

% Package constant physics parameters for the derivative function
derivParams.q_charge = params.q_charge;
derivParams.Ez = params.Ez;
derivParams.E0 = ref_init.E0; % Use REFERENCE E0 for map calculation
derivParams.daOrder = daParams.daOrder;
derivParams.nVars = daParams.nVars;
f_ode = @(z_eval, S_eval) test_state_deriv(z_eval, S_eval, derivParams);

% Integrator options
options = struct();
options.RelTol = 1e-9;       % Tighter tolerance for map calculation
options.AbsTol = 1e-12;
options.MaxStep = rk4Params.dz_step * 10; % Allow larger steps if possible
options.MinStep = rk4Params.dz_step / 1000;
options.InitialStep = rk4Params.dz_step; % Start with previous fixed step
options.Stats = 'on';

% Call the integrator
z_span = [rk4Params.z_start, rk4Params.z_end];
[S_final_map, zf, ~] = integrate_dopri87_da(f_ode, z_span, S_identity, options);

map_calc_time = toc;
fprintf('Finished calculating DA transfer map using DoPri8(7) in %.2f seconds.\n', map_calc_time);
fprintf('Integration reached z = %.6f\n', zf);

% Display the map (optional, can be long)
% fprintf('Final r map:\n');
% S_final_map{1}.displaySeries('%.4e', {'r0', 'pr0'}, 1e-9);
% fprintf('Final pr map:\n');
% S_final_map{2}.displaySeries('%.4e', {'r0', 'pr0'}, 1e-9);

% --- Define Fan of Initial Conditions ---
num_rays_per_dim = 5; % Number of rays = (2*num+1)^2
r0_dev_max = 0.01; % Max deviation in r0 from reference
pr0_dev_max = 0.01; % Max deviation in pr0 from reference

r0_devs = linspace(-r0_dev_max, r0_dev_max, 2*num_rays_per_dim + 1);
pr0_devs = linspace(-pr0_dev_max, pr0_dev_max, 2*num_rays_per_dim + 1);

[R0_dev_grid, PR0_dev_grid] = meshgrid(r0_devs, pr0_devs);
initial_r0s = ref_init.r0 + R0_dev_grid(:);
initial_pr0s = ref_init.pr0 + PR0_dev_grid(:);
num_fan_rays = length(initial_r0s);

fprintf('Evaluating map and analytical solution for %d initial conditions...\n', num_fan_rays);

% --- Store Comparison Results ---
r_final_predicted = zeros(num_fan_rays, 1);
pr_final_predicted = zeros(num_fan_rays, 1);
r_final_analytical = zeros(num_fan_rays, 1);
pr_final_analytical = zeros(num_fan_rays, 1);

% --- Evaluate Map and Analytical Solution for Each Ray ---
az = params.q_charge * params.Ez / params.m_mass; % Constant axial acceleration
z_final = rk4Params.z_end;
z_start = rk4Params.z_start;

for k = 1:num_fan_rays
    r0_k = initial_r0s(k);
    pr0_k = initial_pr0s(k);

    % --- DA Map Prediction ---
    initial_values_k = [r0_k, pr0_k];
    r_final_predicted(k) = evaluate_da_map(S_final_map{1}, initial_values_k);
    pr_final_predicted(k) = evaluate_da_map(S_final_map{2}, initial_values_k);

    % --- Analytical Calculation ---
    % Need vz0_k corresponding to this initial state, assuming same E0 as reference
    % E0 = 0.5*m*(vr0^2 + vz0^2) - q*Ez*z0 => vz0 = sqrt( (2/m)*(E0 + q*Ez*z0) - vr0^2 )
    % where vr0 = pr0 * vz0 => vr0^2 = pr0^2 * vz0^2
    % => vz0^2 = (2/m)*(E0 + q*Ez*z0) - pr0^2 * vz0^2
    % => vz0^2 * (1 + pr0^2) = (2/m)*(E0 + q*Ez*z0)
    vz0_k_sq = (2/params.m_mass)*(ref_init.E0 + params.q_charge*params.Ez*z_start) / (1 + pr0_k^2);

    if vz0_k_sq < 0
        fprintf('Warning: Analytical vz0^2 < 0 for ray %d. Skipping.\n', k);
        r_final_analytical(k) = NaN;
        pr_final_analytical(k) = NaN;
        continue;
    end
    vz0_k = sqrt(vz0_k_sq);
    vr0_k = pr0_k * vz0_k;

    % Calculate time t to reach z_final
    delta_z = z_final - z_start;
    sqrt_term = vz0_k^2 + 2 * az * delta_z;

    if sqrt_term < 0
        fprintf('Warning: Analytical reflection predicted for ray %d. Skipping.\n', k);
        r_final_analytical(k) = NaN;
        pr_final_analytical(k) = NaN;
        continue;
    end

    if abs(az) > eps % Avoid division by zero if no field
        t_final_k = (-vz0_k + sqrt(sqrt_term)) / az;
    else
        t_final_k = delta_z / vz0_k;
    end

    % Calculate final state analytically
    r_final_analytical(k) = r0_k + vr0_k * t_final_k;
    vz_final_k = vz0_k + az * t_final_k;
    if abs(vz_final_k) < eps
        pr_final_analytical(k) = sign(vr0_k) * Inf; % Handle division by zero
    else
        pr_final_analytical(k) = vr0_k / vz_final_k;
    end
end

% --- Plot Comparison ---
figure(3); clf;
subplot(1, 2, 1);
plot(r_final_analytical, r_final_predicted, 'b.');
hold on;
plot(xlim, xlim, 'r--', 'DisplayName', 'Ideal'); % Line y=x
xlabel('Analytical Final r');
ylabel('Predicted Final r (DA Map)');
title('Map Prediction vs Analytical (Final r)');
grid on; axis equal;

subplot(1, 2, 2);
plot(pr_final_analytical, pr_final_predicted, 'b.');
hold on;
plot(xlim, xlim, 'r--', 'DisplayName', 'Ideal'); % Line y=x
xlabel('Analytical Final p_r');
ylabel('Predicted Final p_r (DA Map)');
title('Map Prediction vs Analytical (Final p_r)');
grid on; axis equal;

sgtitle(sprintf('DA Map Validation (Order %d, dz=%.4f)', daParams.daOrder, rk4Params.dz_step));

% --- NEW: Plot Error Heatmaps ---
figure(4); clf;
r_error = r_final_predicted - r_final_analytical;
pr_error = pr_final_predicted - pr_final_analytical;

% Get grid dimensions
n_r_devs = length(r0_devs);
n_pr_devs = length(pr0_devs);

% Reshape errors into grid format
r_error_grid = reshape(r_error, n_pr_devs, n_r_devs); % Note: meshgrid order vs reshape order
pr_error_grid = reshape(pr_error, n_pr_devs, n_r_devs);

% Plot r error heatmap
subplot(1, 2, 1);
imagesc(r0_devs, pr0_devs, r_error_grid);
colorbar;
axis xy; % Place (0,0) deviation at the bottom-left
xlabel('Initial r_0 Deviation');
ylabel('Initial p_{r0} Deviation');
title('Error in Final r (Predicted - Analytical)');
xticks(r0_devs(1:max(1,floor(n_r_devs/4)):n_r_devs)); % Adjust tick frequency if needed
yticks(pr0_devs(1:max(1,floor(n_pr_devs/4)):n_pr_devs));

% Plot pr error heatmap
subplot(1, 2, 2);
imagesc(r0_devs, pr0_devs, pr_error_grid);
colorbar;
axis xy; % Place (0,0) deviation at the bottom-left
xlabel('Initial r_0 Deviation');
ylabel('Initial p_{r0} Deviation');
title('Error in Final p_r (Predicted - Analytical)');
xticks(r0_devs(1:max(1,floor(n_r_devs/4)):n_r_devs));
yticks(pr0_devs(1:max(1,floor(n_pr_devs/4)):n_pr_devs));


sgtitle(sprintf('DA Map Prediction Errors (Order %d, dz=%.4f)', daParams.daOrder, rk4Params.dz_step));

% --- NEW: Plot Overlayed Analytical Trajectories for the Fan ---
figure(5); clf;
hold on;

fprintf('Calculating and plotting analytical trajectories for the fan...\n');

z_traj_points = linspace(z_start, z_final, 200); % Points along z for plotting trajectories
delta_z_traj = z_traj_points - z_start;

colors = cool(num_fan_rays); % Colormap for rays

for k = 1:num_fan_rays
    r0_k = initial_r0s(k);
    pr0_k = initial_pr0s(k);

    % Calculate initial vz0_k and vr0_k consistent with reference E0
    vz0_k_sq = (2/params.m_mass)*(ref_init.E0 + params.q_charge*params.Ez*z_start) / (1 + pr0_k^2);
    if vz0_k_sq < 0, continue; end % Skip if cannot start
    vz0_k = sqrt(vz0_k_sq);
    vr0_k = pr0_k * vz0_k;

    % Calculate analytical trajectory r(z) for this ray
    sqrt_term_traj = vz0_k^2 + 2 * az * delta_z_traj;

    % Find points before potential reflection
    valid_traj_idx = sqrt_term_traj >= 0;
    if ~all(valid_traj_idx)
        first_reflect_idx = find(~valid_traj_idx, 1);
        z_plot = z_traj_points(1:first_reflect_idx-1);
        sqrt_term_plot = sqrt_term_traj(1:first_reflect_idx-1);
    else
        z_plot = z_traj_points;
        sqrt_term_plot = sqrt_term_traj;
    end

    if isempty(z_plot), continue; end % Skip if reflects immediately

    if abs(az) > eps % Avoid division by zero if no field
        t_traj_k = (-vz0_k + sqrt(sqrt_term_plot)) / az;
    else
        t_traj_k = (z_plot - z_start) / vz0_k;
    end

    r_traj_k = r0_k + vr0_k * t_traj_k;

    % Plot this trajectory
    plot(z_plot, r_traj_k, '-', 'Color', [colors(k,:), 0.7], 'LineWidth', 1.0);

end

% Plot the reference trajectory separately for emphasis
ref_sqrt_term_traj = ref_init.vz0^2 + 2 * az * delta_z_traj;
ref_valid_traj_idx = ref_sqrt_term_traj >= 0;
ref_first_reflect_idx = find(~ref_valid_traj_idx, 1);
if isempty(ref_first_reflect_idx), ref_first_reflect_idx = length(z_traj_points)+1; end
ref_z_plot = z_traj_points(1:ref_first_reflect_idx-1);
ref_sqrt_term_plot = ref_sqrt_term_traj(1:ref_first_reflect_idx-1);
if abs(az) > eps
    ref_t_traj = (-ref_init.vz0 + sqrt(ref_sqrt_term_plot)) / az;
else
    ref_t_traj = (ref_z_plot - z_start) / ref_init.vz0;
end
ref_r_traj = ref_init.r0 + ref_init.vr0 * ref_t_traj;
plot(ref_z_plot, ref_r_traj, 'k-', 'LineWidth', 2.0, 'DisplayName', 'Reference Analytical Traj.');


hold off;
xlabel('Axial Position z');
ylabel('Radial Position r');
title('Analytical Trajectories for Initial Condition Fan');
grid on;
axis tight;

% --- End of Updated Plotting Section ---


% ======================================================
% --- Local Helper Functions ---
% ======================================================

function result = evaluate_da_map(da_object, initial_values)
% Evaluates a DA object (power series) for given initial values.
% INPUTS:
%   da_object: A DiffAlg object.
%   initial_values: A row vector [value_var1, value_var2, ...]
%                   Must have length equal to da_object.nVars.
    arguments
        da_object DiffAlg
        initial_values (1,:) double
    end

    if length(initial_values) ~= da_object.nVars
        error('Number of initial values (%d) must match number of DA variables (%d).', ...
              length(initial_values), da_object.nVars);
    end

    coeffs = da_object.Series;
    MI = da_object.getMultiIndices();
    result = 0.0; % Use double precision

    for k = 1:length(coeffs)
        term_coeff = coeffs(k);
        if term_coeff == 0 % Skip zero terms
            continue;
        end

        monomial_value = 1.0; % Value of x1^a1 * x2^a2 * ...
        exponents = MI(k, :);

        for j = 1:da_object.nVars
            if exponents(j) == 0
                continue; % x^0 = 1
            elseif exponents(j) == 1
                monomial_value = monomial_value * initial_values(j);
            else
                monomial_value = monomial_value * (initial_values(j) ^ exponents(j));
            end
        end

        result = result + term_coeff * monomial_value;
    end
end


function dS = test_state_deriv(z_eval, S_eval, params)
% Calculates the state derivative {dr/dz, dp/dz} at z_eval for uniform field.
% (Identical to the function in the previous script)
    % Unpack parameters
    q = params.q_charge;
    Ez = params.Ez;
    E0 = params.E0;
    order = params.daOrder;
    nVars = params.nVars; % Should be 2
    % Unpack state variables (DA objects)
    p = S_eval{2}; % p = dr/dz
    % Calculate scalar Potential Energy U at z_eval
    U_eval = -q * Ez * z_eval - E0;
    % Check for reflection condition U >= 0
    if U_eval >= -1e-12 % Use tolerance
        warning('test_state_deriv:reflection', ...
                'Potential energy U >= 0 (%.3e) at z=%.4f. Particle may reflect.', U_eval, z_eval);
        nan_series = NaN(nchoosek(nVars + order, nVars), 1);
        dS = {DiffAlg(nan_series, order, nVars), DiffAlg(nan_series, order, nVars)};
        return;
    end
    % Calculate scalar dU/dz
    Uz_scalar = -q * Ez;
    % --- Compute Derivatives ---
    one = DiffAlg.one(nVars, order);
    % Eq 1: dr/dz = p
    drdz = p;
    % Eq 2: dp/dz = (1 + p^2) * (dU/dr - p * dU/dz) / (2 * U)
    % Here dU/dr = 0
    one_plus_p2 = one + p * p;
    term_in_parentheses = -1 * (p * Uz_scalar); % Since Ur = 0
    numerator_da = one_plus_p2 * term_in_parentheses;
    % Denominator: 2 * U(z_eval) (scalar)
    denominator_scalar = 2.0 * U_eval;
    if abs(denominator_scalar) < 1e-15
         warning('test_state_deriv:denominator', ...
             'Denominator (2*U) near zero (%.2e) at z=%.3f. Division unstable.', denominator_scalar, z_eval);
         nan_series = NaN(nchoosek(nVars + order, nVars), 1);
         dS = {DiffAlg(nan_series, order, nVars), DiffAlg(nan_series, order, nVars)};
         return;
    end
    % Perform DA division by scalar
    dpdz = numerator_da * (1.0 / denominator_scalar);
    % Return derivatives {dr/dz, dp/dz}
    dS = {drdz, dpdz};
end


% --- RK4 Step Function (identical to the one in run_rk4_da_integration) ---
function S_next = rk4_step_da(z_start_step, h, S_initial, f_deriv)
    % RK4 step where derivative function f_deriv depends on z and S.
    % f_deriv signature: f_deriv(z_eval, S_eval) -> {dr/dz, dp/dz}
    % k1 = h * f(z_start, S_initial)
    try k1_deriv = f_deriv(z_start_step, S_initial);
    catch ME, fprintf('Error in f_deriv (k1) at z=%.4f: %s\n', z_start_step, ME.message); rethrow(ME); end
    if any(isnan(k1_deriv{1}.Series)) || any(isnan(k1_deriv{2}.Series)), S_next = {k1_deriv{1}, k1_deriv{2}}; return; end % Propagate NaN
    k1 = state_scale(k1_deriv, h);
    % k2 = h * f(z_start + h/2, S_initial + k1/2)
    S_temp1 = state_add(S_initial, state_scale(k1, 0.5));
    try k2_deriv = f_deriv(z_start_step + h/2.0, S_temp1);
    catch ME, fprintf('Error in f_deriv (k2) at z=%.4f: %s\n', z_start_step + h/2.0, ME.message); rethrow(ME); end
    if any(isnan(k2_deriv{1}.Series)) || any(isnan(k2_deriv{2}.Series)), S_next = {k2_deriv{1}, k2_deriv{2}}; return; end % Propagate NaN
    k2 = state_scale(k2_deriv, h);
    % k3 = h * f(z_start + h/2, S_initial + k2/2)
    S_temp2 = state_add(S_initial, state_scale(k2, 0.5));
    try k3_deriv = f_deriv(z_start_step + h/2.0, S_temp2);
    catch ME, fprintf('Error in f_deriv (k3) at z=%.4f: %s\n', z_start_step + h/2.0, ME.message); rethrow(ME); end
     if any(isnan(k3_deriv{1}.Series)) || any(isnan(k3_deriv{2}.Series)), S_next = {k3_deriv{1}, k3_deriv{2}}; return; end % Propagate NaN
    k3 = state_scale(k3_deriv, h);
    % k4 = h * f(z_start + h, S_initial + k3)
    S_temp3 = state_add(S_initial, k3);
    try k4_deriv = f_deriv(z_start_step + h, S_temp3);
    catch ME, fprintf('Error in f_deriv (k4) at z=%.4f: %s\n', z_start_step + h, ME.message); rethrow(ME); end
    if any(isnan(k4_deriv{1}.Series)) || any(isnan(k4_deriv{2}.Series)), S_next = {k4_deriv{1}, k4_deriv{2}}; return; end % Propagate NaN
    k4 = state_scale(k4_deriv, h);
    % Combine: S_next = S_initial + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    k_sum = state_add( state_add(k1, state_scale(k2, 2.0)), ...
                       state_add(state_scale(k3, 2.0), k4) );
    S_next = state_add(S_initial, state_scale(k_sum, 1.0/6.0));
end

% --- State Addition Helper ---
function Ssum = state_add(S1, S2)
    Ssum = {S1{1} + S2{1}, S1{2} + S2{2}};
end

% --- State Scaling Helper ---
function Sscaled = state_scale(S, factor)
    Sscaled = {S{1} * factor, S{2} * factor};
end