function [all_R_matrices, z_boundaries, R_total] = run_rk4_da_integration(rk4Params, daParams, bemParams, bemTable)
%run_rk4_da_integration Performs DA integration using RK4 through slices,
%   correctly implementing the trajectory equation including dU/dz.
%
%   INPUTS:
%       rk4Params: Struct with RK4 parameters (z_start, z_end, dz_step)
%       daParams: Struct with DA parameters (daOrder, nVars=2 for state [r, p_r], q_charge, m_mass, E0)
%       bemParams: Struct with BEM results (qVs)
%       bemTable: Table with BEM segment info (r_center, z_center)
%
%   OUTPUTS:
%       all_R_matrices: Cell array of 2x2 transfer matrices for each slice.
%       z_boundaries: Vector of z-coordinates defining slice boundaries.
%       R_total: Overall 2x2 transfer matrix from z_start to z_end.
%
%   Dependencies: DiffAlg class, axial_potential_coeffs function.

    % --- Extract Parameters ---
    z_start = rk4Params.z_start;
    z_end   = rk4Params.z_end;
    dz_step = rk4Params.dz_step;

    daOrder  = daParams.daOrder;
    % Ensure nVars is set correctly for the state vector [r, p_r]
    if daParams.nVars ~= 2
        warning('DA nVars should be 2 for state vector [r, p_r]. Overriding to 2.');
    end
    nVars    = 2; % State variables are r (var 1) and p_r = dr/dz (var 2)

    q_charge = daParams.q_charge;
    m_mass   = daParams.m_mass; % Not used in this dynamics formulation
    E0       = daParams.E0;

    qVs      = bemParams.qVs;
    r_center_bem = bemTable.r_center;
    z_center_bem = bemTable.z_center;

    % --- Slice Definition ---
    total_length = z_end - z_start;
    if total_length <= 0, error('Total z range is non-positive.'); end

    num_full_slices = floor(total_length / dz_step);
    dz_remainder = total_length - num_full_slices * dz_step;
    remainder_tolerance = 1e-9 * dz_step;
    has_remainder = dz_remainder > remainder_tolerance;
    num_total_slices = num_full_slices + has_remainder;

    fprintf('  RK4 Integration Setup (Revamped):\n');
    fprintf('    Range: [%.3f, %.3f], dz_step: %.4f\n', z_start, z_end, dz_step);
    fprintf('    Particle Initial KE (E0): %.4f. Using U = q*V_elec - E0.\n', E0);
    fprintf('    DA Order: %d, nVars (state): %d\n', daOrder, nVars);
    fprintf('    Total slices: %d (%d full + %d remainder).\n', num_total_slices, num_full_slices, has_remainder);

    if num_total_slices == 0, error('No slices defined.'); end

    % --- Initialization ---
    r_da = DiffAlg.var(1, daOrder, nVars); % r (DA var 1)
    pr_da = DiffAlg.var(2, daOrder, nVars);% p_r = dr/dz (DA var 2)

    % Get indices for linear terms extraction
    MI = r_da.getMultiIndices();
    idx_r = find(all(MI == [1, 0], 2), 1); % Index for r^1 * p_r^0
    idx_pr = find(all(MI == [0, 1], 2), 1);% Index for r^0 * p_r^1
    if isempty(idx_r) || isempty(idx_pr), error('Could not find indices for linear terms r or p_r.'); end

    % Store results
    all_R_matrices = cell(num_total_slices, 1);
    z_boundaries = zeros(num_total_slices + 1, 1);
    z_boundaries(1) = z_start;
    R_total = eye(2); % Initialize total transfer matrix

    % Package shared parameters for the derivative function
    derivParams.q_charge = q_charge;
    derivParams.E0 = E0;
    derivParams.daOrder = daOrder;
    derivParams.nVars = nVars; % nVars for state space
    derivParams.bemParams = bemParams; % Contains qVs
    derivParams.bemTable = bemTable; % Contains r_center, z_center_bem

    % --- Loop Through Slices using RK4 ---
    tic;
    current_z_start = z_start; % Start z of the current slice
    S_current = {r_da, pr_da}; % Initial state (identity map)

    for i = 1:num_total_slices
        % Determine current slice thickness
        if i <= num_full_slices
            dz_current_slice = dz_step;
        else % Remainder slice
            dz_current_slice = dz_remainder;
        end
        current_z_end = current_z_start + dz_current_slice;

        % --- Progress Update ---
        if mod(i, 25) == 0 || i == 1 || i == num_total_slices
             fprintf('    Processing slice %d / %d (z: %.3f -> %.3f, dz=%.4f)...\n', ...
                     i, num_total_slices, current_z_start, current_z_end, dz_current_slice);
             if i > 1
                 fprintf('      Traversed z: %.3f\n', current_z_start - z_start);
                 fprintf('      R_total so far: [ % 7.3f  % 7.3f ; % 7.3f  % 7.3f ] (det=%.6f)\n', ...
                         R_total(1,1), R_total(1,2), R_total(2,1), R_total(2,2), det(R_total));
             end
        end

        % --- Define the derivative function handle for this slice ---
        % It needs access to the parameters struct
        f_ode = @(z_eval, S_eval) state_deriv_rk4(z_eval, S_eval, derivParams);

        % --- Perform ONE RK4 step for this slice ---
        % The RK4 step function needs to call f_ode at different z values
        S_next = rk4_step_da(current_z_start, dz_current_slice, S_current, f_ode);

        % --- Extract results for the slice ---
        % The result S_next represents the transformation over the slice.
        % To get the R matrix for the slice, we look at the linear terms
        % of S_next = {r_next_da, p_next_da} relative to the identity {r_da, p_da}.
        r_next_da = S_next{1};
        p_next_da = S_next{2};

        % Extract Slice R Matrix
        R11 = r_next_da.Series(idx_r);  % Coeff of r in r_next
        R12 = r_next_da.Series(idx_pr); % Coeff of p_r in r_next
        R21 = p_next_da.Series(idx_r);  % Coeff of r in p_next
        R22 = p_next_da.Series(idx_pr); % Coeff of p_r in p_next
        R_slice = [R11, R12; R21, R22];
        all_R_matrices{i} = R_slice;

        % Store end z of current slice
        z_boundaries(i+1) = current_z_end;

        % --- Update Total R Matrix ---
        % Note: R_total tracks the map from the *initial* z_start.
        % R_slice transforms from the *start* of the slice to the *end*.
        % R_total(z_end) = R_slice(N) * ... * R_slice(2) * R_slice(1) * R_total(z_start=eye)
        R_total = R_slice * R_total;

        % --- Update state and starting z for the next slice ---
        S_current = S_next; % The output state map becomes the input for the next slice
        current_z_start = current_z_end;

    end % --- End of Slice Loop ---
    loop_time = toc;
    fprintf('    Finished RK4 integration loop (%d slices) in %.2f seconds.\n', num_total_slices, loop_time);

    % --- Display Final Total Matrix ---
    fprintf('  Final Total Transfer Matrix R (z=%.2f to z=%.2f):\n', z_start, z_end);
    fprintf('    [ % 7.3f  % 7.3f ]\n', R_total(1,1), R_total(1,2));
    fprintf('    [ % 7.3f  % 7.3f ]\n', R_total(2,1), R_total(2,2));
    fprintf('    det(R_total) = %.15f\n', det(R_total));

end % --- End of main function ---


% ======================================================
% --- Local Helper Functions for RK4 DA Integration ---
% ======================================================

function dS = state_deriv_rk4(z_eval, S_eval, params)
% Calculates the state derivative {dr/dz, dp/dz} at a specific z_eval.
% Uses DA objects for r and p=dr/dz (nVars=2).
% Constructs necessary potential U(r) and its r-derivative as DA objects
% at z_eval, and calculates dU/dz as a scalar using C1 coefficient at z_eval.

    % Unpack parameters
    q_charge = params.q_charge;
    E0 = params.E0;
    daOrder = params.daOrder;
    nVars = params.nVars; % Should be 2
    qVs = params.bemParams.qVs;
    r_center_bem = params.bemTable.r_center;
    z_center_bem = params.bemTable.z_center;

    % Unpack state variables (DA objects, nVars=2)
    r = S_eval{1};
    p = S_eval{2}; % p = dr/dz

    % 1. Calculate Axial Coefficients for V_elec at z_eval
    Nmax_coeffs = max([daOrder * 2, 10, 1]); % Ensure C1 exists
    Coeffs_eval = axial_potential_coeffs(Nmax_coeffs, z_eval, qVs, r_center_bem, z_center_bem);

    % 2. Construct DA V_elec(r) at z_eval (using only r_da, nVars=2)
    %    Need DA variable 'r' (var 1) for this construction.
    r_da_local = DiffAlg.var(1, daOrder, nVars); % Use state nVars=2
    V_elec_at_z_da = DiffAlg.one(nVars, daOrder) * Coeffs_eval(1); % C0 term
    r_sq_local = r_da_local * r_da_local;
    r_pow_2k_local = r_sq_local;
    max_k_pot = floor(daOrder / 2);

    for k = 1:max_k_pot
        if (2*k + 1) > numel(Coeffs_eval), break; end
        V_2k = factorial(2*k) * Coeffs_eval(2*k + 1); % V_elec^(2k) at z_eval
        coeff_k = ((-1)^k) / ( (4^k) * (factorial(k)^2) );
        term_da = r_pow_2k_local * (coeff_k * V_2k);
        V_elec_at_z_da = V_elec_at_z_da + term_da;
        if k < max_k_pot && (2*(k+1) <= daOrder)
            r_pow_2k_local = r_pow_2k_local * r_sq_local;
        end
    end

    % 3. Construct DA U(r) at z_eval
    U_at_z_da = (V_elec_at_z_da * q_charge) - E0 * DiffAlg.one(nVars, daOrder);

    % 4. Calculate DA dU/dr at z_eval
    %    Differentiate U_at_z_da with respect to r (variable 1)
    Ur_at_z_da = U_at_z_da.differentiate(1);

    % 5. Calculate SCALAR dU/dz at z_eval using C1 coefficient
    if numel(Coeffs_eval) >= 2
        V_elec_z_eval = Coeffs_eval(2); % C1 for V_elec at z_eval
    else
        warning('state_deriv_rk4:coeffs', 'Insufficient C1 coefficient at z=%.3f. Setting Uz to 0.', z_eval);
        V_elec_z_eval = 0;
    end
    Uz_at_z_scalar = q_charge * V_elec_z_eval; % Scalar dU/dz estimate

    % --- Now compute the derivatives using the current state S_eval ---
    one = DiffAlg.one(nVars, daOrder);

    % Eq 1: dr/dz = p
    drdz = p;

    % Eq 2: dp/dz = (1 + p^2) * (dU/dr - p * dU/dz) / (2 * U)
    one_plus_p2 = one + p * p;

    % Evaluate the DA potential U(r) and its derivative dU/dr at the *current state r*
    % This requires substituting r_da_local with the actual state variable r
    U_eval = U_at_z_da;
    Ur_eval = Ur_at_z_da;

    % Check potential energy
    U0_const_eval = U_eval.Series(1);
    if U0_const_eval >= -1e-12
        warning('state_deriv_rk4:potential', 'U >= 0 (%.2e) at z=%.3f. Possible reflection.', U0_const_eval, z_eval);
    end

    % Calculate term in parentheses: dU/dr - p * dU/dz
    % Ur_eval is DA, p is DA, Uz_at_z_scalar is scalar
    term_in_parentheses = Ur_eval - (p * Uz_at_z_scalar);

    % Numerator: (1 + p^2) * (dU/dr - p * dU/dz)
    numerator_da = one_plus_p2 * term_in_parentheses;

    % Denominator: 2 * U(r)
    denominator_da = U_eval * 2.0;

    if abs(denominator_da.Series(1)) < 1e-15
         warning('state_deriv_rk4:denominator', 'Denominator (2*U) near zero (%.2e) at z=%.3f. Division unstable.', denominator_da.Series(1), z_eval);
    end

    % Perform DA division
    dpdz = numerator_da / denominator_da;

    % Return derivatives {dr/dz, dp/dz}
    dS = {drdz, dpdz};

end

% --- RK4 Step Function (passes z_eval to derivative function) ---
function S_next = rk4_step_da(z_start_step, h, S_initial, f_deriv)
    % RK4 step where derivative function f_deriv depends on z and S.
    % f_deriv signature: f_deriv(z_eval, S_eval)

    % k1 = h * f(z_start, S_initial)
    k1_deriv = f_deriv(z_start_step, S_initial);
    k1 = state_scale(k1_deriv, h);

    % k2 = h * f(z_start + h/2, S_initial + k1/2)
    S_temp1 = state_add(S_initial, state_scale(k1, 0.5)); % S_initial + k1/2
    k2_deriv = f_deriv(z_start_step + h/2.0, S_temp1);
    k2 = state_scale(k2_deriv, h);

    % k3 = h * f(z_start + h/2, S_initial + k2/2)
    S_temp2 = state_add(S_initial, state_scale(k2, 0.5)); % S_initial + k2/2
    k3_deriv = f_deriv(z_start_step + h/2.0, S_temp2);
    k3 = state_scale(k3_deriv, h);

    % k4 = h * f(z_start + h, S_initial + k3)
    S_temp3 = state_add(S_initial, k3); % S_initial + k3
    k4_deriv = f_deriv(z_start_step + h, S_temp3);
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
