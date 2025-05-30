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
Vs = [0; 5; 3];
% Charge distribution for a SPECIFIC VOLTAGE SETTING!
qVs = qs' * Vs;

zs = linspace(min(bemTable.z_center), max(bemTable.z_center), 1024);
U_analytical = axial_potential(zs, qVs, bemTable.r_center, bemTable.z_center);
clf;
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

%% Differential Algebra Transfer Map Calculation (Symplectic - Lie Series) - Loop over Slices

% --- DA and Loop Parameters ---
daOrder = 3;  % DA order (adjust as needed, higher order = slower)
daNVars = 2;  % Number of variables (r, p_r)
q_charge = 1.0; % Physical constants
m_mass = 1.0;

% --- NEW: Slice Definition ---
dz_step = 1.0; % Define the desired thickness for the main slices ****** ADJUST THIS ******
z_start = 2;
z_end   = 28;
total_length = z_end - z_start;

if total_length <= 0
    error('Total z range is non-positive. Check zs array.');
end

num_full_slices = floor(total_length / dz_step);
dz_remainder = total_length - num_full_slices * dz_step; % Use subtraction for precision

% Define a tolerance for negligible remainder slice
remainder_tolerance = 1e-9 * dz_step;

has_remainder = dz_remainder > remainder_tolerance;
num_total_slices = num_full_slices + has_remainder;

fprintf('Calculating DA maps using dz_step=%.4f over range [%.2f, %.2f].\n', dz_step, z_start, z_end);
fprintf('Total slices: %d (%d full slices + %d remainder slice).\n', num_total_slices, num_full_slices, has_remainder);

if num_total_slices == 0
     error('No slices defined. Check dz_step and z range.');
end
% --- End NEW Slice Definition ---


% --- Initialization ---
% Initialize DA variables ONCE
r_da = DiffAlg.var(1, daOrder, daNVars);
pr_da = DiffAlg.var(2, daOrder, daNVars);

% Find indices for r and p_r terms ONCE
MI = r_da.getMultiIndices();
idx_r = find(all(MI == [1, 0], 2), 1);  % Index for r^1 * pr^0
idx_pr = find(all(MI == [0, 1], 2), 1); % Index for r^0 * pr^1
if isempty(idx_r) || isempty(idx_pr)
    error('Could not find indices for linear terms (r or p_r) in DA structure.');
end

% Store results for each slice (Optional, but needed for R_total)
all_R_matrices = cell(num_total_slices, 1); % Adjusted size
% all_r_maps = cell(num_total_slices, 1); % Uncomment to store full DA maps
% all_pr_maps = cell(num_total_slices, 1); % Uncomment to store full DA maps

% Constant potential offset
potential_offset = max(abs(U_analytical)) + 1e-6; % Add small buffer

R_total = eye(2); % Initialize total transfer matrix

% --- NEW Loop Through Slices ---
tic; % Start timer
current_z = z_start;
for i = 1:num_total_slices

    % Determine current slice thickness and center
    if i <= num_full_slices
        dz_current_slice = dz_step;
    else % This is the remainder slice
        dz_current_slice = dz_remainder;
    end
    z_center = current_z + dz_current_slice / 2;

    if mod(i, 10) == 0 || i == num_total_slices % Print progress update
         fprintf('  Processing slice %d / %d (z: %.3f -> %.3f, dz=%.4f)...\n', ...
                 i, num_total_slices, current_z, current_z + dz_current_slice, dz_current_slice);
    end

    % 1. Calculate Axial Coefficients for the slice center
    Coeffs_slice = axial_potential_coeffs(Nmax, z_center, qVs, bemTable.r_center, bemTable.z_center);

    % 2. Construct DA Potential V for the slice
    Coeffs_slice(1) = Coeffs_slice(1) + potential_offset; % Apply offset
    if Coeffs_slice(1) <= 0
         warning('Potential + offset is non-positive (%.2e) at z=%.3f. Adjust offset?', Coeffs_slice(1), z_center);
         Coeffs_slice(1) = eps;
    end

    V_da = DiffAlg.one(daNVars, daOrder) * Coeffs_slice(1);
    r_sq = r_da * r_da;
    r_pow_2k = r_sq;
    max_k = floor(min(Nmax, daOrder) / 2);

    for k = 1:max_k
        if (2*k+1) > numel(Coeffs_slice), break; end
        U_2k = factorial(2*k) * Coeffs_slice(2*k + 1);
        coeff_k = ((-1)^k) / ( (4^k) * (factorial(k)^2) );
        V_da = V_da + r_pow_2k * (coeff_k * U_2k);
        if k < max_k && (2*(k+1) <= daOrder)
            r_pow_2k = r_pow_2k * r_sq;
        end
    end

    % 3. Construct DA Hamiltonian H for the slice
    H_da = dsqrt( (pr_da * pr_da) / (2 * m_mass) + V_da * q_charge );

    % 4. Calculate Lie Series Map for the slice step dz_current_slice
    g_da = H_da * (-dz_current_slice); % Use current slice thickness
    r_new_da_slice = r_da;
    pr_new_da_slice = pr_da;
    r_term = r_da;
    pr_term = pr_da;

    for k_lie = 1:daOrder
        r_term_next_pb = poisson_bracket(r_term, g_da);
        pr_term_next_pb = poisson_bracket(pr_term, g_da);
        r_term_k = r_term_next_pb * (1/k_lie);
        pr_term_k = pr_term_next_pb * (1/k_lie);
        r_new_da_slice = r_new_da_slice + r_term_k;
        pr_new_da_slice = pr_new_da_slice + pr_term_k;
        r_term = r_term_k;
        pr_term = pr_term_k;
    end

    % 5. Extract Slice R Matrix
    R11 = r_new_da_slice.Series(idx_r);
    R12 = r_new_da_slice.Series(idx_pr);
    R21 = pr_new_da_slice.Series(idx_r);
    R22 = pr_new_da_slice.Series(idx_pr);
    R_slice = [R11, R12; R21, R22];
    all_R_matrices{i} = R_slice;

    % 6. Update Total R Matrix
    R_total = R_slice * R_total;

    % 7. Update starting z for the next slice
    current_z = current_z + dz_current_slice;

end % --- End of NEW Slice Loop ---
loop_time = toc;
fprintf('Finished calculating maps for %d slices in %.2f seconds.\n', num_total_slices, loop_time);

% --- Display Final Results ---
fprintf('\nTotal First-Order Transfer Matrix R (Product of all slices):\n');
fprintf('[ % 7.3f  % 7.3f ]\n', R_total(1,1), R_total(1,2));
fprintf('[ % 7.3f  % 7.3f ]\n', R_total(2,1), R_total(2,2));
fprintf('\ndet(R_total) = %.15f\n', det(R_total));