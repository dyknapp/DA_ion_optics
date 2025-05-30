% Parameters
res = 256;
spacing = 1;
z0 = 0;    % starting z position
R = 1;

% Define a cell array of electrode definitions.
% Each definition is a structure with:
%   - a shape function (with a function handle that takes res as input)
%   - an optional thickness field (used to update z position for the next electrode)
electrodeData = {
    struct('thickness', 10, 'shape', @(res) tube(R, 10, res))
    struct('thickness', 10, 'shape', @(res) tube(R, 10, res))
    struct('thickness', 10, 'shape', @(res) tube(R, 10, res))
};

% Build the electrode vector automatically.
electrodes = createElectrodeStack(electrodeData, 'z0', z0, 'spacing', spacing, 'res', res);

% BEM analysis
[qs, bemTable] = BEM_monopole(electrodes);

%% Voltage configuration
Vs = [0; -0.5; -1.0];
% Charge distribution for a SPECIFIC VOLTAGE SETTING!
qVs = qs' * Vs;

zs = linspace(min(bemTable.z_center), max(bemTable.z_center), 1024);
U_analytical = axial_potential(zs, qVs, bemTable.r_center, bemTable.z_center);
figure(1); clf;
plot(zs, U_analytical, '-k')

Nmax = 10; % Choose the maximum order for the series expansion
z0_expansion = 12.1;
expansion_range = 0.5;

% Ensure qVs is a column vector (if needed)
if ~iscolumn(qVs); qVs = qVs.'; end

% Calculate the coefficients C_0, C_1, ..., C_Nmax
Coeffs = axial_potential_coeffs(Nmax, z0_expansion, qVs, bemTable.r_center, bemTable.z_center);

% Display the coefficients
disp('Axial Potential Coefficients C_n:');
for n = 0:Nmax
    fprintf('C_%3d = % .6e\n', n, Coeffs(n+1)); % Coeffs(n+1) corresponds to C_n
end

local_idcs = (zs >= z0_expansion - expansion_range) & (zs <= z0_expansion + expansion_range);
z_local = zs(local_idcs);
delta_z = z_local - z0_expansion;
U_approx = zeros(size(z_local));
for n = 0:Nmax
    U_approx = U_approx + Coeffs(n+1) * delta_z.^n;
end
hold on
plot(z_local, U_approx, '-r');
xline(z0_expansion, '-.r');
hold off
title('Approximated Potential from Series'); xlabel('z'); ylabel('U(z)')

fprintf('Mean abs rel error: %.6e\n', ...
    mean(abs(U_approx - U_analytical(local_idcs)) ./ abs(U_analytical(local_idcs))))

Cs = zeros([3 numel(zs)]);
for idx = 1:numel(zs)
    Cs(:, idx) = axial_potential_coeffs(2, zs(idx), qVs, bemTable.r_center, bemTable.z_center);
end
hold on
yyaxis right
plot(zs, Cs(3, :), '-b')
hold off

legend('Analytical BEM potential', 'Series expansion', 'Expansion Center', 'Refractive power')

% Differential Algebra Transfer Map Calculation (Symplectic - Lie Series)
% Define DA parameters (same as before)
daOrder = 10;  % Example DA order, adjust as needed
daNVars = 2;  % Number of variables (r, p_r)

% Define physical constants (same as before)
q_charge = 1.0;
m_mass = 1.0;

% Initialize DA variables (same as before)
r_da = DiffAlg.var(1, daOrder, daNVars);  % DA variable for r
pr_da = DiffAlg.var(2, daOrder, daNVars); % DA variable for p_r

% Construct the potential V(r, z0_expansion) as a DA object (same as before)
% The potential should be shifted by a constant so that it's non-negative.
% The easiest approach is just to set the constant term to a big number.
Coeffs(1) = max(max(abs(U_analytical)), eps);
V_da = DiffAlg.one(daNVars, daOrder) * Coeffs(1);
r_sq = r_da * r_da; % Precompute r^2
r_pow_2k = r_sq;    % Initialize r^(2k) for k=1
max_k = floor(Nmax / 2); % Use Nmax from the script
for k = 1:max_k
    if 2*k > daOrder, break; end
    U_2k = factorial(2*k) * Coeffs(2*k + 1);
    coeff_k = ((-1)^k) / ( (4^k) * (factorial(k)^2) );
    V_da = V_da + r_pow_2k * (coeff_k * U_2k); %
    if k < max_k
        r_pow_2k = r_pow_2k * r_sq; %
    end
end

% Construct the Hamiltonian H = sqrt(pr^2/(2m) + qV) as a DA object (same as before)
H_da = dsqrt( (pr_da * pr_da) / (2 * m_mass) + V_da * q_charge ); % Using dsqrt from DiffAlg

% --- Symplectic Integration using Lie Series ---
% Map M = exp(-dz * L_H), where L_H f = {f, H}
% Apply map: (r_new, pr_new) = M * (r, pr) = exp(L_g) * (r, pr)
% where g = -dz * H

dz = expansion_range; % Step size for this example slice
g_da = H_da * (-dz);  % The generator for the step

% Initialize result with the identity map (0th order term)
r_new_da = r_da;      %
pr_new_da = pr_da;    %

% Initialize the current term in the series (k=0 term)
r_term = r_da;        %
pr_term = pr_da;      %


function pb = poisson_bracket(f, g)
    % Computes the Poisson bracket {f, g} for two DiffAlg objects f, g
    % Assumes nVars = 2, where variable 1 is r and variable 2 is p_r.
    arguments
        f DiffAlg
        g DiffAlg
    end
    if f.Order ~= g.Order || f.nVars ~= g.nVars
        error('Poisson bracket requires DA objects of the same Order and nVars.');
    end
    if f.nVars ~= 2
        error('Poisson bracket as defined here requires nVars = 2 (r, p_r).');
    end

    % Calculate partial derivatives
    df_dr = f.differentiate(1);   % df/dr
    dg_dpr = g.differentiate(2);  % dg/dp_r
    df_dpr = f.differentiate(2);  % df/dp_r
    dg_dr = g.differentiate(1);   % dg/dr

    % Compute the bracket: (df/dr * dg/dpr) - (df/dpr * dg/dr)
    pb = (df_dr * dg_dpr) - (df_dpr * dg_dr); %
end

% Calculate Lie series: f_new = sum_{k=0}^{order} {f, g}_k / k!
% where {f, g}_0 = f, {f, g}_{k+1} = {{f, g}_k, g}
fprintf('Calculating Lie series terms (Order %d):\n', daOrder);
for k = 1:daOrder
    fprintf('  Term k =%3d: calculating Poisson brackets...\n', k);
    % Calculate the next term: {current_term, g}
    r_term_next_pb = poisson_bracket(r_term, g_da); %
    pr_term_next_pb = poisson_bracket(pr_term, g_da); %

    % Divide by k for the series term {f, g}_k / k! = ({f, g}_{k-1}/(k-1)!, g) / k
    r_term_k = r_term_next_pb * (1/k); %
    pr_term_k = pr_term_next_pb * (1/k); %

    % Add the term to the sum
    r_new_da = r_new_da + r_term_k; %
    pr_new_da = pr_new_da + pr_term_k; %

    % Update the term for the next iteration
    r_term = r_term_k;
    pr_term = pr_term_k;

    % Optional: Check if terms become negligible (e.g., norm of series coeffs)
    % if norm(r_term.Series) < 1e-15 && norm(pr_term.Series) < 1e-15
    %     fprintf('  Terms negligible, stopping series early at k=%d.\n', k);
    %     break;
    % end
end
fprintf('Lie series calculation complete.\n');

% Display the resulting transfer map (optional)
fprintf('Symplectic Transfer Map for slice z0=%.2f, dz=%.2f (Order %d):\n', z0_expansion, dz, daOrder);
fprintf('\nNew r(z0+dz):');
r_new_da.displaySeries('% .6f', {'r0', 'p_r,0'}, 1.0e-6);
fprintf('\nNew p_r(z0+dz):');
pr_new_da.displaySeries('% .6f', {'r0', 'p_r,0'}, 1.0e-6);

% --- Extract First-Order Transfer Matrix ---

fprintf('\nExtracting First-Order Transfer Matrix (R)...\n');

% Ensure we are working with 2 variables (r, p_r)
if r_new_da.nVars ~= 2 || pr_new_da.nVars ~= 2
    error('Transfer matrix extraction requires nVars = 2.');
end

% Get the multi-index mapping (should be the same for r_new_da and pr_new_da)
MI = r_new_da.getMultiIndices();

% Find the index corresponding to the 'r' term (multi-index [1, 0])
idx_r = find(all(MI == [1, 0], 2), 1);

% Find the index corresponding to the 'p_r' term (multi-index [0, 1])
idx_pr = find(all(MI == [0, 1], 2), 1);

% Check if linear terms were found (they should exist if Order >= 1)
if isempty(idx_r) || isempty(idx_pr)
    error('Could not find indices for linear terms (r or p_r) in DA map.');
end

% Extract the coefficients from the Series vectors
R11 = r_new_da.Series(idx_r);
R12 = r_new_da.Series(idx_pr);
R21 = pr_new_da.Series(idx_r);
R22 = pr_new_da.Series(idx_pr);

% Assemble the matrix
R_matrix = [R11, R12; R21, R22];

% Print the matrix
fprintf('First-Order Transfer Matrix R:\n');
fprintf('[ % .9e  % .9e ]\n', R_matrix(1,1), R_matrix(1,2));
fprintf('[ % .9e  % .9e ]\n', R_matrix(2,1), R_matrix(2,2));

fprintf('\ndet(M) = %.15f\n', det(R_matrix))

%% Differential Algebra Transfer Map Calculation (RK4 Integration) - Loop over Slices

% --- DA and Loop Parameters ---
daOrder = 10;
nVars = 2;
q_charge = 1.0; % Set charge (+1.0 for proton, -1.0 for electron)
m_mass = 1.0;
E0 = 1.0; % Initial Kinetic Energy (where V_elec = 0) ****** ADJUST THIS ******

% --- Slice Definition ---
dz_step = 0.01; % ****** ADJUST THIS ******
z_start = 3;
z_end   = 28;
total_length = z_end - z_start;
xline([z_start, z_end], '--k')

if total_length <= 0, error('Total z range is non-positive.'); end

num_full_slices = floor(total_length / dz_step);
dz_remainder = total_length - num_full_slices * dz_step;
remainder_tolerance = 1e-9 * dz_step;
has_remainder = dz_remainder > remainder_tolerance;
num_total_slices = num_full_slices + has_remainder;

fprintf('Calculating DA maps using RK4 with dz_step=%.4f over range [%.2f, %.2f].\n', dz_step, z_start, z_end);
fprintf('Particle Initial KE (E0) = %.4f. Using U = q*V_elec - E0.\n', E0);
fprintf('Total slices: %d (%d full slices + %d remainder slice).\n', num_total_slices, num_full_slices, has_remainder);

if num_total_slices == 0, error('No slices defined.'); end

% --- Initialization ---
r_da = DiffAlg.var(1, daOrder, nVars);
pr_da = DiffAlg.var(2, daOrder, nVars);
MI = r_da.getMultiIndices();
idx_r = find(all(MI == [1, 0], 2), 1);
idx_pr = find(all(MI == [0, 1], 2), 1);
if isempty(idx_r) || isempty(idx_pr), error('Could not find indices for linear terms.'); end
all_R_matrices = cell(num_total_slices, 1);
R_total = eye(2);

% --- Loop Through Slices using RK4 ---
tic;
current_z = z_start;
for i = 1:num_total_slices
    % Determine current slice thickness and center
    if i <= num_full_slices
        dz_current_slice = dz_step;
    else % Remainder slice
        dz_current_slice = dz_remainder;
    end
    z_center = current_z + dz_current_slice / 2;

    if mod(i, 25) == 0 || i == num_total_slices % Print progress update less frequently
         fprintf('  Processing slice %d / %d (z: %.3f -> %.3f, dz=%.4f)...\n', ...
                 i, num_total_slices, current_z, current_z + dz_current_slice, dz_current_slice);
         fprintf('    Traversed z-distance so far: % .3f\n', current_z - z_start)
         fprintf('    Ray transfer matrix  so far:\n')
         fprintf('    [ % 7.3f  % 7.3f ]\n', R_total(1,1), R_total(1,2));
         fprintf('    [ % 7.3f  % 7.3f ]\n', R_total(2,1), R_total(2,2));
         fprintf('    det(R_total) = %.15f\n', det(R_total));
         fprintf('\n');
    end

    % 1. Calculate Axial Coefficients for V_elec
    Nmax_coeffs = max(Nmax, daOrder*2);
    Coeffs_slice = axial_potential_coeffs(Nmax_coeffs, z_center, qVs, bemTable.r_center, bemTable.z_center);

    % 2. Construct DA Electrostatic Potential V_elec for the slice
    V_elec_slice_da = DiffAlg.one(daNVars, daOrder) * Coeffs_slice(1); % V_elec C0
    r_sq = r_da * r_da;
    r_pow_2k = r_sq;
    max_k_pot = floor(min(Nmax_coeffs, daOrder) / 2);

    for k = 1:max_k_pot
        if (2*k+1) > numel(Coeffs_slice), break; end
        V_2k = factorial(2*k) * Coeffs_slice(2*k + 1); % V_elec^(2k)
        coeff_k = ((-1)^k) / ( (4^k) * (factorial(k)^2) );
        term_da = r_pow_2k * (coeff_k * V_2k);
        V_elec_slice_da = V_elec_slice_da + term_da;
        if k < max_k_pot && (2*(k+1) <= daOrder)
            r_pow_2k = r_pow_2k * r_sq;
        end
    end

    % 3. Construct U = q*V_elec - E0 for the dynamics equation
    %    **** MODIFIED based on E0 convention ****
    U_slice_da = (V_elec_slice_da * q_charge) - E0 * DiffAlg.one(nVars, daOrder);

    % *** Check if U (potential energy, should be negative) is non-negative ***
    % This corresponds to KE = -U becoming non-positive (particle stops/reflects)
    if U_slice_da.Series(1) >= -1e-12 % Check constant term: q*V_elec(z_center, r=0) - E0
        warning('Potential energy U = qV-E0 is non-negative (%.2e) at z=%.3f, r=0. Particle may reflect.', U_slice_da.Series(1), z_center);
        % Depending on desired behavior, you might stop integration,
        % flag the result, or let the numerics proceed (might error).
        % For now, we proceed but issue a warning.
    end
    % *** Check if U is near zero for division ***
    if abs(U_slice_da.Series(1)) < 1e-12
        warning('Potential energy U=qV-E0 is near zero (%.2e) at z=%.3f. Division in dp/dz might fail.', U_slice_da.Series(1), z_center);
    end

    % 4. Define the state derivative function handle
    f_slice = @(S) state_deriv_adapted(S, U_slice_da); % Pass U = qV-E0

    % 5. Initialize state for the slice integration (Identity Map)
    S_initial_slice = {r_da, pr_da};

    % 6. Perform ONE RK4 step for this slice
    S_final_slice = rk4_step(S_initial_slice, dz_current_slice, f_slice);

    % 7. Extract results for the slice
    r_new_da_slice = S_final_slice{1};
    pr_new_da_slice = S_final_slice{2};

    % 8. Extract Slice R Matrix
    R11 = r_new_da_slice.Series(idx_r);
    R12 = r_new_da_slice.Series(idx_pr);
    R21 = pr_new_da_slice.Series(idx_r);
    R22 = pr_new_da_slice.Series(idx_pr);
    R_slice = [R11, R12; R21, R22];
    all_R_matrices{i} = R_slice;

    % --- Add code to store z boundaries ---
    if i == 1
        z_boundaries = zeros(num_total_slices + 1, 1);
        z_boundaries(1) = z_start;
    end
    z_boundaries(i+1) = current_z + dz_current_slice; % Store end z of current slice

    % 9. Update Total R Matrix
    R_total = R_slice * R_total;

    % 10. Update starting z for the next slice
    current_z = current_z + dz_current_slice;

end % --- End of Slice Loop ---
loop_time = toc;
fprintf('Finished calculating maps for %d slices using RK4 in %.2f seconds.\n', num_total_slices, loop_time);

% --- Display Final Results ---
fprintf('\nTotal First-Order Transfer Matrix R (Product of all slices):\n');
fprintf('[ % 7.3f  % 7.3f ]\n', R_total(1,1), R_total(1,2));
fprintf('[ % 7.3f  % 7.3f ]\n', R_total(2,1), R_total(2,2));
fprintf('\ndet(R_total) = %.15f\n', det(R_total));

% --- Ray fan plot ---

fprintf('\nPerforming ray tracing for visualization...\n');

% --- Define Initial Ray Fan ---
num_rays = 11; % Number of rays to trace (odd number includes central ray)
max_r0 = R * 0.5; % Maximum initial radius (e.g., 50% of tube radius R)
initial_r0s = linspace(-max_r0, max_r0, num_rays);
initial_p0 = 0; % Initial slope (dr/dz)
initial_states = [initial_r0s; repmat(initial_p0, 1, num_rays)]; % 2 x num_rays

% --- Trace Rays Through Slices ---
num_steps = num_total_slices;
r_positions = zeros(num_steps + 1, num_rays); % Store r at each boundary
p_positions = zeros(num_steps + 1, num_rays); % Store p at each boundary (optional)

% Set initial conditions
r_positions(1, :) = initial_states(1, :);
p_positions(1, :) = initial_states(2, :);
current_states = initial_states;

% Apply R matrix for each slice
for i = 1:num_steps
    R_slice = all_R_matrices{i};
    current_states = R_slice * current_states; % Apply map: X_new = R * X_old
    r_positions(i+1, :) = current_states(1, :);
    p_positions(i+1, :) = current_states(2, :);
end

% --- Plot Ray Fan ---
figure(2); clf; % Create a new figure window
hold on;

% Plot the traced rays
plot(z_boundaries, r_positions, '.-b'); % Plot r vs z for all rays

% Plot electrode geometry for context
for k = 1:numel(electrodes)
    elec = electrodes(k);
    plot(elec.zs, elec.rs, '-k', 'LineWidth', 1.5);
    plot(elec.zs, -elec.rs, '-k', 'LineWidth', 1.5);
end

% Add labels and title
xlabel('Axial Position z');
ylabel('Radial Position r');
title(sprintf('Ray Fan Plot (Linear Tracing, %d Rays, p0=%.2f)', num_rays, initial_p0));
grid on;
axis tight; % Adjust axis limits to data
ylim_curr = ylim; % Get current y-limits
ylim_max = max(abs(ylim_curr)); % Find max absolute y
ylim([-max(ylim_max, R * 1.1), max(ylim_max, R * 1.1)]); % Set symmetric y-limits, at least +/- tube radius
hold off;

fprintf('Ray fan plot generated.\n');

%% === Parameters for Mesh Generation ===
R_outer = 50; % Outer radius of the surrounding electrode blocks (adjust as needed)
num_angular_steps = 64; % Controls rotational smoothness
cutaway_angle_degrees = 90; % Angle to remove (90 for 1/4 cutaway = 3/4 remaining)
output_folder = 'electrode_and_traj_meshes_stl'; % Folder to save STL files
combined_stl_filename = sprintf('./%s/combined.stl', output_folder);
export_trajectory = true;

% Convert cutaway angle to radians and define the angular range
cutaway_angle_rad = deg2rad(cutaway_angle_degrees);
theta = linspace(0, 2*pi - cutaway_angle_rad, num_angular_steps);
if abs(theta(end) - 2*pi) < 1e-9 && cutaway_angle_degrees == 0
    theta = theta(1:end-1);
    num_angular_steps = num_angular_steps -1;
end

% Ensure output folder exists
if ~exist(output_folder, 'dir')
   mkdir(output_folder);
end

% === Initialization for Combined Geometry ===
all_vertices = []; % Stores vertices from all components
all_faces = [];    % Stores faces from all components (with offset indices)
vertex_offset = 0; % Running count of vertices added so far

% === Part 1: Loop Through Electrodes to Generate Blocks ===
num_electrodes = numel(electrodes);
fprintf('\n--- Processing Electrode Blocks for Combined File ---\n');
for i = 1:num_electrodes
    fprintf('Processing electrode block %d...\n', i);

    % --- Get Profile Points ---
    profile_z = electrodes(i).zs;
    r_inner = electrodes(i).rs;
    profile_z = profile_z(:);
    r_inner = r_inner(:);

    if numel(profile_z) ~= numel(r_inner)
        warning('Skipping electrode block %d: Mismatch between number of z points (%d) and r points (%d).', i, numel(profile_z), numel(r_inner));
        continue;
    end
    num_profile_points = numel(profile_z);
    if num_profile_points < 2
        warning('Skipping electrode block %d: Needs at least 2 profile points (found %d).', i, num_profile_points);
        continue;
    end

    % Define outer radius profile
    r_outer = ones(num_profile_points, 1) * R_outer;

    % --- Generate 3D Vertices for Current Block ---
    N_inner_points_per_angle = num_profile_points;
    N_total_inner_vertices = N_inner_points_per_angle * num_angular_steps;
    N_total_outer_vertices = N_inner_points_per_angle * num_angular_steps;
    current_block_vertices = zeros(N_total_inner_vertices + N_total_outer_vertices, 3); % Preallocate for current block only
    all_inner_vertices_current = zeros(N_total_inner_vertices, 3); % Temporary for calculation
    all_outer_vertices_current = zeros(N_total_outer_vertices, 3); % Temporary for calculation

    for k = 1:num_angular_steps
        angle_idx_start_in_array = (k-1) * N_inner_points_per_angle;
        indices_for_this_angle = angle_idx_start_in_array + (1:N_inner_points_per_angle);
        current_angle = theta(k);
        cos_theta = cos(current_angle);
        sin_theta = sin(current_angle);
        x_inner = r_inner .* cos_theta;
        y_inner = r_inner .* sin_theta;
        all_inner_vertices_current(indices_for_this_angle, :) = [x_inner, y_inner, profile_z];
        x_outer = r_outer .* cos_theta;
        y_outer = r_outer .* sin_theta;
        all_outer_vertices_current(indices_for_this_angle, :) = [x_outer, y_outer, profile_z];
    end
    current_block_vertices = [all_inner_vertices_current; all_outer_vertices_current];
    outer_vertices_start_idx_local = N_total_inner_vertices + 1; % Local start index within current_block_vertices
    N_vertices_current = size(current_block_vertices, 1);

    % --- Generate Faces for Current Block (Indices relative to current_block_vertices) ---
    faces_current_block = []; % Initialize faces for THIS block
    faces_inner_cyl = []; faces_outer_cyl = []; faces_end_cap1 = []; faces_end_cap2 = []; faces_cut_surf1 = []; faces_cut_surf2 = []; % Clear sub-arrays

    for k = 1:(num_angular_steps - 1)
        for j = 1:(num_profile_points - 1)
            idx1_in = (k-1)*num_profile_points + j;
            idx2_in = k*num_profile_points + j;
            idx3_in = k*num_profile_points + (j+1);
            idx4_in = (k-1)*num_profile_points + (j+1);
            faces_inner_cyl = [faces_inner_cyl; idx1_in, idx3_in, idx2_in; idx1_in, idx4_in, idx3_in];

            idx1_out = outer_vertices_start_idx_local + (idx1_in - 1); % Use local start index
            idx2_out = outer_vertices_start_idx_local + (idx2_in - 1);
            idx3_out = outer_vertices_start_idx_local + (idx3_in - 1);
            idx4_out = outer_vertices_start_idx_local + (idx4_in - 1);
            faces_outer_cyl = [faces_outer_cyl; idx1_out, idx2_out, idx3_out; idx1_out, idx3_out, idx4_out];
        end
        j_start = 1; j_end = num_profile_points;
        idx1_in_cap1 = (k-1)*num_profile_points + j_start;
        idx2_in_cap1 = k*num_profile_points + j_start;
        idx1_out_cap1 = outer_vertices_start_idx_local + (idx1_in_cap1 - 1);
        idx2_out_cap1 = outer_vertices_start_idx_local + (idx2_in_cap1 - 1);
        faces_end_cap1 = [faces_end_cap1; idx1_in_cap1, idx2_out_cap1, idx2_in_cap1; idx1_in_cap1, idx1_out_cap1, idx2_out_cap1];

        idx4_in_cap2 = (k-1)*num_profile_points + j_end;
        idx3_in_cap2 = k*num_profile_points + j_end;
        idx4_out_cap2 = outer_vertices_start_idx_local + (idx4_in_cap2 - 1);
        idx3_out_cap2 = outer_vertices_start_idx_local + (idx3_in_cap2 - 1);
        faces_end_cap2 = [faces_end_cap2; idx4_in_cap2, idx3_in_cap2, idx3_out_cap2; idx4_in_cap2, idx3_out_cap2, idx4_out_cap2];
    end
    k_start = 1; k_end = num_angular_steps;
    for j = 1:(num_profile_points - 1)
        idx1_in_cut1 = (k_start-1)*num_profile_points + j;
        idx4_in_cut1 = (k_start-1)*num_profile_points + (j+1);
        idx1_out_cut1 = outer_vertices_start_idx_local + (idx1_in_cut1 - 1);
        idx4_out_cut1 = outer_vertices_start_idx_local + (idx4_in_cut1 - 1);
        faces_cut_surf1 = [faces_cut_surf1; idx1_in_cut1, idx1_out_cut1, idx4_out_cut1; idx1_in_cut1, idx4_out_cut1, idx4_in_cut1];

        idx2_in_cut2 = (k_end-1)*num_profile_points + j;
        idx3_in_cut2 = (k_end-1)*num_profile_points + (j+1);
        idx2_out_cut2 = outer_vertices_start_idx_local + (idx2_in_cut2 - 1);
        idx3_out_cut2 = outer_vertices_start_idx_local + (idx3_in_cut2 - 1);
        faces_cut_surf2 = [faces_cut_surf2; idx2_in_cut2, idx3_in_cut2, idx3_out_cut2; idx2_in_cut2, idx3_out_cut2, idx2_out_cut2];
    end
    faces_current_block = [faces_inner_cyl; faces_outer_cyl; faces_end_cap1; faces_end_cap2; faces_cut_surf1; faces_cut_surf2];

    % --- Append to Global Lists ---
    if ~isempty(faces_current_block) && ~isempty(current_block_vertices)
        fprintf('  Appending %d vertices and %d faces for block %d.\n', N_vertices_current, size(faces_current_block, 1), i);
        all_vertices = [all_vertices; current_block_vertices]; % Append vertices
        all_faces = [all_faces; faces_current_block + vertex_offset]; % Append faces OFFSET by current vertex count
        vertex_offset = vertex_offset + N_vertices_current; % UPDATE offset for the NEXT block/trajectory
    else
        warning('  Skipping append for block %d due to empty vertices or faces.', i);
    end

end % End of electrode block loop
fprintf('--- Finished Electrode Block Processing ---\n');
fprintf('Total vertices after electrode blocks: %d\n', vertex_offset);


% === Part 2: Generate and Append Ion Trajectory Solid ===
z_traj = z_boundaries;
r_traj = r_positions(:, end);
if export_trajectory && ~isempty(z_traj) % Ensure z_traj is valid if export_trajectory is true
    fprintf('\n--- Processing Ion Trajectory Solid for Combined File ---\n');
    num_traj_profile_points = numel(z_traj);
    if num_traj_profile_points < 2
        warning('Skipping trajectory solid: Needs at least 2 trajectory points.');
    else
        % --- Generate 3D Vertices for Trajectory ---
        traj_vertices_current = zeros(num_traj_profile_points * num_angular_steps, 3);
        center_vertices_current = [zeros(num_traj_profile_points, 1), zeros(num_traj_profile_points, 1), z_traj];
        num_traj_surf_verts = size(traj_vertices_current, 1);
        center_vertices_start_idx_local = num_traj_surf_verts + 1; % Local index within vertices_traj_combined

        for k = 1:num_angular_steps
             angle_idx_start = (k-1)*num_traj_profile_points;
             indices_for_this_angle = angle_idx_start + (1:num_traj_profile_points);
             current_angle = theta(k);
             x_traj = r_traj .* cos(current_angle);
             y_traj = r_traj .* sin(current_angle);
             traj_vertices_current(indices_for_this_angle, :) = [x_traj, y_traj, z_traj];
        end
        vertices_traj_combined = [traj_vertices_current; center_vertices_current]; % Combine surface and center points for THIS trajectory part
        N_vertices_traj = size(vertices_traj_combined, 1);

        % --- Generate Faces for Trajectory (Indices relative to vertices_traj_combined) ---
        faces_traj_current = []; % Initialize faces for THIS trajectory
        faces_traj_curved = []; faces_traj_cap1 = []; faces_traj_cap2 = []; faces_traj_cut1 = []; faces_traj_cut2 = []; % Clear sub-arrays

        for k = 1:(num_angular_steps - 1)
            fprintf('%3d/%3d\n', k, (num_angular_steps - 1))
            for j = 1:(num_traj_profile_points - 1)
                idx1 = (k-1)*num_traj_profile_points + j;
                idx2 = k*num_traj_profile_points + j;
                idx3 = k*num_traj_profile_points + (j+1);
                idx4 = (k-1)*num_traj_profile_points + (j+1);
                faces_traj_curved = [faces_traj_curved; idx1, idx2, idx3; idx1, idx3, idx4];
            end
            j_start = 1; j_end = num_traj_profile_points;
            center_idx_start = center_vertices_start_idx_local + j_start -1;
            center_idx_end = center_vertices_start_idx_local + j_end -1;
            idx_outer_1_cap1 = (k-1)*num_traj_profile_points + j_start;
            idx_outer_2_cap1 = k*num_traj_profile_points + j_start;
            faces_traj_cap1 = [faces_traj_cap1; idx_outer_1_cap1, center_idx_start, idx_outer_2_cap1];
            idx_outer_1_cap2 = (k-1)*num_traj_profile_points + j_end;
            idx_outer_2_cap2 = k*num_traj_profile_points + j_end;
            faces_traj_cap2 = [faces_traj_cap2; idx_outer_1_cap2, idx_outer_2_cap2, center_idx_end];
        end
         k_start = 1; k_end = num_angular_steps;
         for j = 1:(num_traj_profile_points - 1)
             idx_outer_1_cut1 = (k_start-1)*num_traj_profile_points + j;
             idx_outer_2_cut1 = (k_start-1)*num_traj_profile_points + (j+1);
             idx_inner_1_cut1 = center_vertices_start_idx_local + j - 1;
             idx_inner_2_cut1 = center_vertices_start_idx_local + j;
             faces_traj_cut1 = [faces_traj_cut1; idx_inner_1_cut1, idx_outer_1_cut1, idx_outer_2_cut1; idx_inner_1_cut1, idx_outer_2_cut1, idx_inner_2_cut1];

             idx_outer_1_cut2 = (k_end-1)*num_traj_profile_points + j;
             idx_outer_2_cut2 = (k_end-1)*num_traj_profile_points + (j+1);
             idx_inner_1_cut2 = center_vertices_start_idx_local + j - 1;
             idx_inner_2_cut2 = center_vertices_start_idx_local + j;
             faces_traj_cut2 = [faces_traj_cut2; idx_inner_1_cut2, idx_outer_2_cut2, idx_outer_1_cut2; idx_inner_1_cut2, idx_inner_2_cut2, idx_outer_2_cut2];
         end
        faces_traj_current = [faces_traj_curved; faces_traj_cap1; faces_traj_cap2; faces_traj_cut1; faces_traj_cut2];

        % --- Append to Global Lists ---
        if ~isempty(faces_traj_current) && ~isempty(vertices_traj_combined)
            fprintf('  Appending %d vertices and %d faces for trajectory.\n', N_vertices_traj, size(faces_traj_current, 1));
            all_vertices = [all_vertices; vertices_traj_combined]; % Append vertices
            all_faces = [all_faces; faces_traj_current + vertex_offset]; % Append faces OFFSET by current vertex count
            vertex_offset = vertex_offset + N_vertices_traj; % Update offset (though not used after this)
        else
             warning('  Skipping append for trajectory due to empty vertices or faces.');
        end
    end % End if num_traj_profile_points check
    fprintf('--- Finished Ion Trajectory Processing ---\n');
else
     fprintf('\n--- Skipped Ion Trajectory Processing ---\n');
end % End if export_trajectory check

% === Part 3: Final Export of Combined Geometry ===
fprintf('\n--- Exporting Combined Geometry --- \n');
fprintf('Total combined vertices: %d\n', size(all_vertices, 1));
fprintf('Total combined faces: %d\n', size(all_faces, 1));

if isempty(all_vertices) || isempty(all_faces)
    warning('No geometry generated. Skipping final STL export.');
else
    % --- Create Combined Triangulation Object ---
    fprintf('Creating final triangulation object...\n');
    try
        T_combined = triangulation(all_faces, all_vertices);
        fprintf('Triangulation object created successfully.\n');

        % --- Export Combined Geometry to STL ---
        fprintf('Writing combined STL file: %s...\n', combined_stl_filename);
        try
            stlwrite(T_combined, combined_stl_filename);
            fprintf('Successfully wrote combined STL file.\n');
        catch ME_write
             warning('Failed to write combined STL file %s. Error: %s', combined_stl_filename, ME_write.message);
        end % End try-catch stlwrite

    catch ME_tri
        warning('Failed to create final combined triangulation object. Error: %s', ME_tri.message);
        fprintf('    Check for issues like non-manifold geometry or inconsistent face definitions.\n');
    end % End try-catch triangulation
end % End if isempty check

disp('--- Script Finished ---');


%%
% ==================================================
% --- Helper Functions (Need to be defined below or accessible) ---
% ==================================================

function dS = state_deriv_adapted(S, U_slice)
    % Adapted state derivative for dr/dz = p, dp/dz = (1+p^2)*U_r / (2 U)
    % U_slice is a DA object representing the potential U(r) for the current slice.
    % It depends only on r (variable 1).
    r = S{1}; p = S{2};
    one = DiffAlg.one(r.nVars, r.Order); % Assumes DiffAlg class is accessible

    % dr/dz = p
    drdz = p;

    % dp/dz = (1+p^2)*U_r / (2 U)
    p2 = one + p * p;
    U_r = U_slice.differentiate(1); % Derivative w.r.t r (var 1)
    % U_z is zero because U_slice is constructed from Coeffs at fixed z_center

    % Ensure denominator is not zero (handled by potential_offset usually)
    if abs(U_slice.Series(1)) < 1e-12 % Check constant term magnitude
         warning('Constant term of U_slice is near zero (%.2e) during derivative calculation.', U_slice.Series(1));
         % Consider returning zero derivative or erroring if this happens
    end

    % Use division overload (mrdivide) or inverse: dpdz = (p2 * U_r) * inverse(U_slice * 2)
     dpdz = (p2 * U_r) / (U_slice * 2.0);

    dS = {drdz, dpdz};
end

% --- RK4 Step Function ---
function S_next = rk4_step(S, h, f)
    % Standard 4th-order Runge-Kutta step for state S = {r, p}
    % S is cell array {DA_obj, DA_obj}
    % h is step size
    % f is function handle f(S) returning {dr/dz, dp/dz} as DA objects
    k1 = f(S); % {drdz1, dpdz1}
    S_temp = state_add(S, state_scale(k1, h/2));
    k2 = f(S_temp); % {drdz2, dpdz2}
    S_temp = state_add(S, state_scale(k2, h/2));
    k3 = f(S_temp); % {drdz3, dpdz3}
    S_temp = state_add(S, state_scale(k3, h));
    k4 = f(S_temp); % {drdz4, dpdz4}

    % Combine steps: S + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    k_sum = state_add( state_add(k1, state_scale(k2, 2.0)), ...
                       state_add(state_scale(k3, 2.0), k4) );
    S_next = state_add(S, state_scale(k_sum, h/6.0));
end

% --- State Addition Helper ---
function Ssum = state_add(S1, S2)
    % Adds two states S1={r1,p1}, S2={r2,p2} elementwise
    Ssum = {S1{1} + S2{1}, S1{2} + S2{2}};
end

% --- State Scaling Helper ---
function Sscaled = state_scale(S, factor)
    % Multiplies state S={r,p} by a scalar factor
    Sscaled = {S{1} * factor, S{2} * factor};
end
