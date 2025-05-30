% Script to test and visualize the Gaussian_potential_expansion function,
% initialize DA objects, compute Hamiltonian K and Poisson brackets,
% calculate Lie transfer maps for multiple slices, perform linear ray tracing,
% and include visualization plots.

clear; clc; close all;
set(0,'DefaultFigureWindowStyle','docked') % Dock figures for better management

% --- Parameters ---
z0_expansion = 0.0;      % Point for initial DA potential construction (used for plots)
Nmax_expansion = 10;     % Maximum order for expansions
daOrder = Nmax_expansion;% DA order matches expansion order
z_plot_half_range = 0.5; % Plot range for z in potential plots
r_plot_max = 0.5;        % Maximum r value for radial/3D plots
num_points_plot = 500;   % Number of points for 2D plotting
num_ribs_3d = 15;        % Number of radial "ribs" to plot in 3D
num_points_r_rib = 100;  % Number of points for each radial rib plot

% --- Integration and Tracing Parameters ---
z_start_trace = -10.0;
z_end_trace = 10.0;
dz_trace = 0.025;
z_boundaries_trace = z_start_trace:dz_trace:z_end_trace;
num_slices_trace = numel(z_boundaries_trace) - 1;

% --- Ray Fan Parameters ---
traceParams.num_rays = 15; % Number of rays (odd recommended)
traceParams.max_r0 = 0.1;  % Max initial radius
traceParams.initial_p0 = 0.0; % Initial p_r = m*dr/dt

fprintf('Gaussian Potential Simulation using Lie Maps (Order %d)\n', daOrder);
fprintf('Integration Range: z = [%.2f, %.2f], dz = %.3f (%d slices)\n', ...
        z_start_trace, z_end_trace, dz_trace, num_slices_trace);
fprintf('Ray Fan: %d rays, max_r0 = %.3f, p_r0 = %.3f\n', ...
        traceParams.num_rays, traceParams.max_r0, traceParams.initial_p0);

% --- Calculate Coefficients (at z0, primarily for plotting) ---
fprintf('\nCalculating coefficients at z0 = %.2f (for plots)...\n', z0_expansion);
try
    [z_coeffs_plot, ~] = Gaussian_potential_expansion(z0_expansion, Nmax_expansion);
    fprintf('Successfully calculated axial coefficients for plots.\n');
catch ME
    fprintf('Error calling Gaussian_potential_expansion: %s\n', ME.message);
    return;
end

% --- Precompute Factorials ---
fact_vals = zeros(Nmax_expansion + 1, 1);
if Nmax_expansion >= 0, fact_vals(1) = 1; end % 0!
for n = 1:Nmax_expansion, fact_vals(n+1) = fact_vals(n) * n; if ~isfinite(fact_vals(n+1)), warning('FactorialOverflow', 'Factorial overflow occurred at n=%d.', n); end; end

% --- DA Initialization (Variables only) ---
fprintf('\nInitializing DA variables...\n');
% Define physical constants
m_mass = 1.0; q_charge = 1.0; E_total = 2.0; % Adjusted E_total for potentially more interesting dynamics
% Define DA parameters
nVars = 3; r_idx = 1; z_idx = 2; pr_idx = 3; varNames_da = {'r', 'z', 'p_r'};
% Initialize canonical DA variables
try
    r_da = DiffAlg.var(r_idx, daOrder, nVars); z_da = DiffAlg.var(z_idx, daOrder, nVars); pr_da = DiffAlg.var(pr_idx, daOrder, nVars); one_da = DiffAlg.one(nVars, daOrder);
    dim_da = nchoosek(nVars + daOrder, nVars); zero_da = DiffAlg(zeros(dim_da, 1), daOrder, nVars);
    % Get indices for linear term extraction
    MI = r_da.getMultiIndices();
    idx_r_lin = find(all(MI == [1, 0, 0], 2), 1); % r^1 z^0 pr^0
    idx_pr_lin = find(all(MI == [0, 0, 1], 2), 1); % r^0 z^0 pr^1
    if isempty(idx_r_lin) || isempty(idx_pr_lin), error('Could not find indices for linear terms r or p_r.'); end
catch ME_da_init, fprintf('Error initializing DA variables: %s\n', ME_da_init.message); return; end
fprintf('DA variables initialized (Order=%d, nVars=%d).\n', daOrder, nVars);


% --- Multi-Slice Lie Map Integration & Linear Matrix Extraction ---
fprintf('\n--- Calculating Transfer Matrix for Each Slice ---\n');

all_R_matrices_trace = cell(num_slices_trace, 1);
R_total_trace = eye(2); % Accumulates the product of linear matrices

integration_successful = true;
tic; % Start timer for integration loop

for i = 1:num_slices_trace
    z_slice_start = z_boundaries_trace(i);
    z_slice_end = z_boundaries_trace(i+1);
    z_center_slice = (z_slice_start + z_slice_end) / 2.0;

    if mod(i, 5) == 0 || i == 1 || i == num_slices_trace
        fprintf('  Processing slice %d/%d (z: %.3f -> %.3f, center: %.3f)...\n', ...
                i, num_slices_trace, z_slice_start, z_slice_end, z_center_slice);
    end

    % --- Recalculate DA potential and Hamiltonian based on z_center_slice ---
    % 1. Coefficients at z_center
    try [z_coeffs_slice, ~] = Gaussian_potential_expansion(z_center_slice, Nmax_expansion);
    catch ME_coeff_slice, fprintf('ERROR slice %d coeff: %s\n', i, ME_coeff_slice.message); integration_successful = false; break; end

    % 2. V(z) DA object around z_center
    V_axial_da_slice = zero_da; z_center_da = one_da * z_center_slice; delta_z_da_slice = z_da - z_center_da; delta_z_pow_n = one_da;
    for n = 0:Nmax_expansion, coeff_n = z_coeffs_slice(n + 1); if isnan(coeff_n) || coeff_n == 0, continue; end; if n > 0, delta_z_pow_n = delta_z_pow_n * delta_z_da_slice; end; V_axial_da_slice = V_axial_da_slice + (delta_z_pow_n * coeff_n); end

    % 3. V^(m)(z) DA objects
    V_derivs_da_slice = cell(Nmax_expansion + 1, 1); V_deriv_da_current = V_axial_da_slice; V_derivs_da_slice{1} = V_deriv_da_current; diff_ok_slice = true;
    for m = 1:Nmax_expansion
        try V_deriv_da_current = V_deriv_da_current.differentiate(z_idx); V_derivs_da_slice{m + 1} = V_deriv_da_current;
        catch ME_diff_slice, fprintf('ERROR slice %d diff m=%d: %s\n', i, m, ME_diff_slice.message); diff_ok_slice = false; break; end
    end; if ~diff_ok_slice, integration_successful = false; break; end

    % 4. V(r, z) DA object
    V_rz_da_slice = zero_da; max_k_rz = floor(Nmax_expansion / 2); r_pow_2k_da = one_da; constr_ok_slice = true; r_sq_da_local = [];
    for k = 0:max_k_rz
        idx_V2k = 2*k + 1; if idx_V2k > numel(V_derivs_da_slice), constr_ok_slice = false; break; end
        V_2k_da = V_derivs_da_slice{idx_V2k}; F_k = calculate_Fk(k, fact_vals); if isnan(F_k) || any(isnan(V_2k_da.Series)), continue; end; if F_k == 0, continue; end
        if k > 0, if k == 1, r_sq_da_local = r_da * r_da; r_pow_2k_da = r_sq_da_local; else, r_pow_2k_da = r_pow_2k_da * r_sq_da_local; end; end
        try term_k_da = (r_pow_2k_da * F_k) * V_2k_da; V_rz_da_slice = V_rz_da_slice + term_k_da;
        catch ME_term_slice, fprintf('ERROR slice %d Vrz term k=%d: %s\n', i, k, ME_term_slice.message); constr_ok_slice = false; break; end
    end; if ~constr_ok_slice, integration_successful = false; break; end

    % 5. Term_DA for slice
    Term_DA_slice = (V_rz_da_slice * (-q_charge) + (one_da * E_total)) * (2 * m_mass);

    % 6. Hamiltonian K for slice
    hamiltonian_ok_slice = false; K_da_slice = zero_da;
    
    try
        pr_sq_da = pr_da * pr_da; InsideSqrt_DA_slice = Term_DA_slice - pr_sq_da; const_term_inside_sqrt = InsideSqrt_DA_slice.Series(1);
        % --- DEBUG: Compare constant term ---
        expected_const = 2*m_mass*(E_total - q_charge*exp(-z_center_slice^2));
        actual_const = InsideSqrt_DA_slice.Series(1);
        if abs(expected_const - actual_const) > 1e-6 % Or some other tolerance
            fprintf('DEBUG Slice %d (z=%.3f): Constant term mismatch! Expected %.6e, Got %.6e\n', i, z_center_slice, expected_const, actual_const);
        end
        % --- End DEBUG ---
        if const_term_inside_sqrt <= 1e-12, warning('Slice %d: Constant term non-positive (%.3e). Using R=eye.', i, const_term_inside_sqrt); K_da_slice = zero_da * NaN; % Mark as invalid
        else pz_da_slice = InsideSqrt_DA_slice.dsqrt(); K_da_slice = pz_da_slice * (-1.0); hamiltonian_ok_slice = true; end
    catch ME_ham_slice, fprintf('ERROR constructing K slice %d: %s\n', i, ME_ham_slice.message); end
    if ~hamiltonian_ok_slice, K_da_slice = zero_da * NaN; end % Ensure invalid if error or non-positive

    % 7. Generator g for slice step dz_trace
    g_da_slice = K_da_slice * (-dz_trace); % Will be NaN if K_da_slice is NaN

    % 8. Lie map for slice
    r_map_slice = zero_da; pr_map_slice = zero_da; map_calc_ok = false;
    if hamiltonian_ok_slice % Only calculate map if K was ok
        try
            r_map_slice = calculate_lie_transfer_map(r_da, g_da_slice, daOrder, fact_vals);
            pr_map_slice = calculate_lie_transfer_map(pr_da, g_da_slice, daOrder, fact_vals);
            map_calc_ok = true; % Check for NaNs below
            if any(isnan(r_map_slice.Series)) || any(isnan(pr_map_slice.Series))
                 warning('Slice %d: NaN detected in calculated map.', i);
                 map_calc_ok = false;
            end
        catch ME_map_slice
            fprintf('ERROR calculating map for slice %d: %s\n', i, ME_map_slice.message);
        end
    end

    % 9. Extract Linear Matrix R_slice
    if map_calc_ok
        R11 = r_map_slice.Series(idx_r_lin); R12 = r_map_slice.Series(idx_pr_lin);
        R21 = pr_map_slice.Series(idx_r_lin); R22 = pr_map_slice.Series(idx_pr_lin);
        R_slice = [R11, R12; R21, R22];
        % Check determinant
        detR = det(R_slice);
        if abs(detR - 1.0) > 1e-3 % Tolerance for numerical errors
            warning('Slice %d: Determinant of R is %.6f (deviates from 1).', i, detR);
            % Optional: Normalize R_slice? R_slice = R_slice / sqrt(abs(detR));
        end
    else
        fprintf('  Slice %d: Using identity matrix due to previous errors.\n', i);
        R_slice = eye(2); % Use identity if map calculation failed
    end

    % Store R_slice and update R_total
    all_R_matrices_trace{i} = R_slice;
    R_total_trace = R_slice * R_total_trace; % R_total(z_end) = R_N * ... * R_1

end % End of slice loop

integration_time = toc;

% Proceed only if loop finished without breaking
if i == num_slices_trace || ~integration_successful % Check if loop completed fully
    fprintf('\nFinished integration loop (%d slices processed) in %.2f seconds.\n', i, integration_time);
    if integration_successful
         fprintf('Total Linear Transfer Matrix R_total (Product order R_N*...*R_1):\n');
         fprintf('  [ % 9.6f  % 9.6f ]\n', R_total_trace(1,1), R_total_trace(1,2));
         fprintf('  [ % 9.6f  % 9.6f ]\n', R_total_trace(2,1), R_total_trace(2,2));
         fprintf('  det(R_total) = %.9f\n', det(R_total_trace));
    else
         fprintf('Integration loop did not complete successfully. Total matrix may be incorrect.\n');
    end

    % --- Linear Ray Tracing ---
    fprintf('\nPerforming linear ray tracing...\n');

    initial_r0s = linspace(-traceParams.max_r0, traceParams.max_r0, traceParams.num_rays);
    initial_states_trace = [initial_r0s; repmat(traceParams.initial_p0, 1, traceParams.num_rays)]; % [r0; pr0]

    num_steps_trace = num_slices_trace;
    r_positions_trace = zeros(num_steps_trace + 1, traceParams.num_rays);
    pr_positions_trace = zeros(num_steps_trace + 1, traceParams.num_rays); % For potential phase space plot

    r_positions_trace(1, :) = initial_states_trace(1, :);
    pr_positions_trace(1, :) = initial_states_trace(2, :);
    current_states_trace = initial_states_trace;

    for slice_idx = 1:num_steps_trace
        R_slice = all_R_matrices_trace{slice_idx};
        if any(isnan(R_slice(:))) % Check if matrix is valid
             warning('NaN detected in R_slice matrix for slice %d. Stopping trace.', slice_idx);
             r_positions_trace(slice_idx+1:end, :) = NaN; % Mark rest as invalid
             pr_positions_trace(slice_idx+1:end, :) = NaN;
             break; % Stop tracing
        end
        current_states_trace = R_slice * current_states_trace; % Apply map: X_out = R_slice * X_in
        r_positions_trace(slice_idx+1, :) = current_states_trace(1, :);
        pr_positions_trace(slice_idx+1, :) = current_states_trace(2, :);
    end

    % --- Plot Ray Fan ---
    fig_ray = figure(3); % Use figure 3 for ray fan
    fig_ray.Name = 'Ray Fan (Lie Map - Linear Trace)';
    plot(z_boundaries_trace, r_positions_trace, 'b-');
    hold on;
    plot(z_boundaries_trace, zeros(size(z_boundaries_trace)), 'k--', 'LineWidth', 0.5); % Axis line
    xlabel('Axial Position (z)');
    ylabel('Radial Position (r)');
    title(sprintf('Ray Fan (%d rays, r0=%.3f, pr0=%.3f, E=%.1f)', traceParams.num_rays, traceParams.max_r0, traceParams.initial_p0, E_total));
    grid on;
    axis tight;
    ylim_curr = ylim;
    ylim_max_abs = max(abs(ylim_curr));
    ylim([-max(ylim_max_abs*1.1, traceParams.max_r0*1.5), max(ylim_max_abs*1.1, traceParams.max_r0*1.5)]); % Symmetric y-axis, ensure initial fan visible
    hold off;
    fprintf('Ray fan plot generated in Figure 3.\n');

else
    fprintf('\nIntegration loop failed or was interrupted. Skipping ray tracing and plotting.\n');
end


% --- 2D Plotting Section (Axial and Radial at z0) ---
% [ Code remains the same ]
fprintf('\nGenerating 2D potential plots...\n');
fig2d = figure(1); fig2d.Name = sprintf('2D Gaussian Potential Expansion (z0=%.1f, Nmax=%d)', z0_expansion, Nmax_expansion); z_values_2d = linspace(z0_expansion - z_plot_half_range, z0_expansion + z_plot_half_range, num_points_plot); V_exact_z_2d = exp(-z_values_2d.^2); V_approx_z_2d = zeros(size(z_values_2d)); delta_z_2d = z_values_2d - z0_expansion; for n = 0:Nmax_expansion, V_approx_z_2d = V_approx_z_2d + z_coeffs_plot(n + 1) * (delta_z_2d.^n); end; rel_error_z_2d = abs(V_approx_z_2d - V_exact_z_2d) ./ abs(V_exact_z_2d); rel_error_z_2d(abs(V_exact_z_2d) < eps) = abs(V_approx_z_2d(abs(V_exact_z_2d) < eps)); r_values_2d = linspace(0, r_plot_max, num_points_plot); V_approx_r_at_z0 = zeros(size(r_values_2d)); max_k_2d = floor(Nmax_expansion / 2); V_derivs_at_z0 = calculate_V_derivatives(z0_expansion, z0_expansion, z_coeffs_plot, Nmax_expansion, fact_vals); for k = 0:max_k_2d; V_2k_at_z0 = V_derivs_at_z0(2*k + 1); F_k = calculate_Fk(k, fact_vals); if isnan(F_k) || isnan(V_2k_at_z0), V_approx_r_at_z0(:) = NaN; break; end; V_approx_r_at_z0 = V_approx_r_at_z0 + F_k * V_2k_at_z0 * (r_values_2d.^(2*k)); end; subplot(1, 2, 1); plot(z_values_2d, V_exact_z_2d, 'k-', 'LineWidth', 2, 'DisplayName', 'Exact V(z) = exp(-z^2)'); hold on; plot(z_values_2d, V_approx_z_2d, 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('Taylor Approx. (Nmax=%d)', Nmax_expansion)); xline(z0_expansion, 'b-.', 'LineWidth', 1, 'DisplayName', sprintf('Expansion Point z0=%.1f', z0_expansion)); hold off; xlabel('Axial Position (z)'); ylabel('Potential V(z, r=0)'); title('Axial Potential: Exact vs. Taylor Series'); legend('show', 'Location', 'best'); grid on; ylim_current = ylim; ylim([min(ylim_current(1), -0.1), max(ylim_current(2), 1.1)]); yyaxis right; semilogy(z_values_2d, rel_error_z_2d, 'm:', 'LineWidth', 1, 'DisplayName', 'Relative Error'); ylabel('Relative Error |Approx - Exact| / |Exact|'); ax = gca; ax.YAxis(2).Color = 'm'; ylim_err = ylim; ylim([max(1e-16, ylim_err(1)), min(1e2, ylim_err(2))]); legend('show', 'Location', 'best'); subplot(1, 2, 2); plot(r_values_2d, V_approx_r_at_z0, 'b-', 'LineWidth', 2); xlabel('Radial Position (r)'); ylabel(sprintf('Potential V(r, z=%.1f)', z0_expansion)); title(sprintf('Radial Potential Approx. at z=%.1f (Nmax=%d)', z0_expansion, Nmax_expansion)); grid on; xlim([0, r_plot_max]); fprintf('2D plotting complete.\n');

% --- 3D Plotting Section ("Boat Hull") ---
% [ Code remains the same ]
fprintf('\nGenerating 3D potential plot...\n'); fig3d = figure(2); fig3d.Name = sprintf('3D Potential V(r,z) Approx. (z0=%.1f, Nmax=%d)', z0_expansion, Nmax_expansion); hold on; z_values_3d_keel = linspace(z0_expansion - z_plot_half_range, z0_expansion + z_plot_half_range, 100); z_rib_values = linspace(min(z_values_3d_keel), max(z_values_3d_keel), num_ribs_3d); r_values_rib = linspace(-r_plot_max, r_plot_max, num_points_r_rib); V_keel_approx = zeros(size(z_values_3d_keel)); delta_z_keel = z_values_3d_keel - z0_expansion; for n = 0:Nmax_expansion, V_keel_approx = V_keel_approx + z_coeffs_plot(n + 1) * (delta_z_keel.^n); end; plot3(zeros(size(z_values_3d_keel)), z_values_3d_keel, V_keel_approx, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Axial Potential (Keel)'); max_k_3d = floor(Nmax_expansion / 2); for i = 1:num_ribs_3d; z_rib = z_rib_values(i); V_derivs_at_rib = calculate_V_derivatives(z_rib, z0_expansion, z_coeffs_plot, Nmax_expansion, fact_vals); V_rib = zeros(size(r_values_rib)); valid_rib = true; for k = 0:max_k_3d; V_2k_at_rib = V_derivs_at_rib(2*k + 1); F_k = calculate_Fk(k, fact_vals); if isnan(F_k) || isnan(V_2k_at_rib), warning('NaN encountered calculating rib at z=%.2f, k=%d. Skipping rib.', z_rib, k); valid_rib = false; break; end; r_pow_2k = r_values_rib.^(2*k); if k==0, r_pow_2k(r_values_rib==0) = 1; end; V_rib = V_rib + F_k * V_2k_at_rib * r_pow_2k; end; if valid_rib, plot3(r_values_rib, ones(size(r_values_rib)) * z_rib, V_rib, 'b-', 'LineWidth', 1, 'HandleVisibility','off'); end; end; hold off; grid on; xlabel('Radial Position (r)'); ylabel('Axial Position (z)'); zlabel('Potential V(r,z)'); title(sprintf('3D Potential Approx. V(r,z) (Nmax=%d)', Nmax_expansion)); view(30, 25); axis tight; legend('show'); fprintf('3D plotting complete.\n');


% =========================================================================
%                            LOCAL FUNCTIONS
% =========================================================================
function [z_coeffs, r_coeffs] = Gaussian_potential_expansion(z0, Nmax) % ... (full function code) ...
    arguments; z0 (1,1) double {mustBeFinite}; Nmax (1,1) double {mustBeInteger, mustBeNonnegative}; end; Y = zeros(Nmax + 1, 1); Y(1) = exp(-z0^2); if Nmax >= 1, Y(2) = -2 * z0 * Y(1); end; for n = 1:(Nmax-1); if ~isfinite(Y(n+1)) || ~isfinite(Y(n)), Y(n + 2:end) = NaN; break; end; Y(n + 2) = -2 * z0 * Y(n + 1) - 2 * n * Y(n); end; z_coeffs = zeros(Nmax + 1, 1); fact_vals = zeros(Nmax + 1, 1); if Nmax >= 0, fact_vals(1) = 1; end; for n = 1:Nmax, fact_vals(n+1) = fact_vals(n) * n; end; for n = 0:Nmax; if fact_vals(n + 1) == 0, z_coeffs(n + 1) = NaN; elseif ~isfinite(Y(n + 1)), z_coeffs(n + 1) = NaN; else, z_coeffs(n + 1) = Y(n + 1) / fact_vals(n + 1); end; end; max_k = floor(Nmax / 2); r_coeffs = zeros(max_k + 1, 1); for k = 0:max_k; idx_C2k = 2*k + 1; C_2k = z_coeffs(idx_C2k); if ~isfinite(C_2k), r_coeffs(k + 1) = NaN; continue; end; F_k = calculate_Fk(k, fact_vals); if isnan(F_k), r_coeffs(k + 1) = NaN; continue; end; r_coeffs(k + 1) = F_k * C_2k; end;
end
function V_derivs = calculate_V_derivatives(z_eval, z0, C_coeffs, Nmax, fact_vals) % ... (full function code) ...
    V_derivs = zeros(Nmax + 1, 1); delta_z = z_eval - z0; for m = 0:Nmax; sum_val = 0; for n = m:Nmax; C_n = C_coeffs(n + 1); if isnan(C_n) || C_n == 0, continue; end; n_perm_m = 1; for i = 0:(m-1), n_perm_m = n_perm_m * (n - i); end; pow_delta_z = delta_z^(n - m); term = C_n * n_perm_m * pow_delta_z; if ~isfinite(term), sum_val = NaN; break; end; sum_val = sum_val + term; end; V_derivs(m + 1) = sum_val; end;
end
function F_k = calculate_Fk(k, fact_vals) % ... (full function code) ...
    idx_2k = 2*k + 1; idx_k = k + 1; if idx_2k > numel(fact_vals) || idx_k > numel(fact_vals), F_k = NaN; return; end; fact_2k_val = fact_vals(idx_2k); fact_k_val = fact_vals(idx_k); pow_4k = 4^k; denom_Fk = pow_4k * (fact_k_val^2); if denom_Fk == 0, F_k = NaN; if ((-1)^k * fact_2k_val) ~= 0, F_k = Inf * sign((-1)^k * fact_2k_val); end; elseif ~isfinite(fact_2k_val) || ~isfinite(denom_Fk), F_k = NaN; else, F_k = ((-1)^k * fact_2k_val) / denom_Fk; end;
end
function pb = poisson_bracket(f, g) % ... (full function code) ...
    arguments; f DiffAlg; g DiffAlg; end; if ~isa(f, 'DiffAlg') || ~isa(g, 'DiffAlg'), error('Inputs must be DiffAlg objects.'); end; if f.Order ~= g.Order || f.nVars ~= g.nVars, error('Poisson bracket requires DA objects of the same Order and nVars.'); end; if f.nVars ~= 3, error('Poisson bracket as defined here requires nVars = 3 (r=1, z=2, p_r=3).'); end; r_idx = 1; pr_idx = 3; try; df_dr  = f.differentiate( r_idx); dg_dpr = g.differentiate(pr_idx); df_dpr = f.differentiate(pr_idx); dg_dr  = g.differentiate( r_idx); catch ME_diff; error('Error during differentiation for Poisson bracket: %s', ME_diff.message); end; pb = (df_dr * dg_dpr) - (df_dpr * dg_dr);
end
function f_final = calculate_lie_transfer_map(f_initial, g, order, fact_vals) % ... (full function code) ...
    arguments; f_initial DiffAlg; g DiffAlg; order (1,1) double {mustBeInteger, mustBeNonnegative}; fact_vals (:,1) double; end; if ~isa(f_initial, 'DiffAlg') || ~isa(g, 'DiffAlg'), error('Inputs must be DiffAlg objects.'); end; if f_initial.Order ~= g.Order || f_initial.nVars ~= g.nVars || order > f_initial.Order, error('Incompatible DA objects or order for Lie map calculation.'); end; if numel(fact_vals) < (order + 1), error('Precomputed factorials vector is too short.'); end; if any(isnan(g.Series)), warning('LieMap:NaN_Generator', 'Generator g contains NaN values. Result will be NaN.'); nan_series = NaN(size(f_initial.Series)); f_final = DiffAlg(nan_series, f_initial.Order, f_initial.nVars); return; end; f_final = f_initial; current_pb = f_initial; for k = 1:order; try; if any(isnan(current_pb.Series)), warning('LieMap:NaN_Propagation', 'NaN detected in intermediate Poisson bracket at k=%d. Stopping series.', k-1); if ~any(isnan(f_final.Series)), nan_series = NaN(size(f_initial.Series)); f_final = DiffAlg(nan_series, f_initial.Order, f_initial.nVars); end; return; end; current_pb = poisson_bracket(current_pb, g); catch ME_pb_map, fprintf('Error calculating Poisson bracket inside Lie map at k=%d: %s\n', k, ME_pb_map.message); nan_series = NaN(size(f_initial.Series)); f_final = DiffAlg(nan_series, f_initial.Order, f_initial.nVars); return; end; fact_k = fact_vals(k + 1); if fact_k == 0 || ~isfinite(fact_k), warning('LieMap:FactorialIssue', 'Factorial(%d) is zero or non-finite (%.2e). Map inaccurate.', k, fact_k); term_k = current_pb * NaN; else, term_k = current_pb * (1.0 / fact_k); end; f_final = f_final + term_k; term_norm = norm(term_k.Series); if ~isfinite(term_norm), warning('LieMap:NonFiniteTerm', 'Term norm is non-finite at k=%d. Stopping series.', k); if ~any(isnan(f_final.Series)), nan_series = NaN(size(f_initial.Series)); f_final = DiffAlg(nan_series, f_initial.Order, f_initial.nVars); end; return; end; if term_norm < 1e-15, fprintf('  Lie series converged early at k=%d (Term norm: %.2e)\n', k, term_norm); break; end; end;
end
