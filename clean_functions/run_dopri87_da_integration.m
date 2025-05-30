function [z_steps, S_steps, stats] = run_dopri87_da_integration(dopriParams, daParams, bemParams, bemTable, integratorOptions)
%run_dopri87_da_integration Performs DA integration using adaptive Dormand-Prince 8(7)
%   over the full specified range [z_start, z_end], returning intermediate steps.
%
%   INPUTS:
%       dopriParams: Struct with integration range (z_start, z_end).
%       daParams: Struct with DA parameters (daOrder, nVars=2 for state [r, p_r], q_charge, m_mass, E0)
%       bemParams: Struct with BEM results (qVs)
%       bemTable: Table with BEM segment info (r_center, z_center)
%       integratorOptions: (Optional) Struct with integrator options passed to
%                          integrate_dopri87_da (e.g., RelTol, AbsTol, Stats).
%
%   OUTPUTS:
%       z_steps   : Vector [z_start, z_step1, ..., z_final] of z values at each accepted step.
%       S_steps   : Cell array {S_initial, S_step1, ..., S_final} of DA state maps at each z_step.
%       stats     : Struct with integration statistics (if requested in options).
%
%   Dependencies: DiffAlg class, axial_potential_coeffs function,
%                 integrate_dopri87_da function, get_dopri87_coeffs function.

    % --- Extract Parameters ---
    z_start = dopriParams.z_start;
    z_end   = dopriParams.z_end;
    z_span = [z_start, z_end]; % Integration span

    daOrder  = daParams.daOrder;
    if daParams.nVars ~= 2
        warning('DA nVars should be 2 for state vector [r, p_r]. Overriding to 2.');
    end
    nVars    = 2; % State variables are r (var 1) and p_r = dr/dz (var 2)

    q_charge = daParams.q_charge;
    E0       = daParams.E0;

    % --- Default Integrator Options ---
    defaultOpts = struct();
    defaultOpts.RelTol = 1e-7;
    defaultOpts.AbsTol = 1e-10;
    defaultOpts.Stats = 'off';
    verbose = true;
    if isfield(dopriParams, 'InitialStep'), defaultOpts.InitialStep = dopriParams.InitialStep; end
    if isfield(dopriParams, 'MaxStep'), defaultOpts.MaxStep = dopriParams.MaxStep; end
    if isfield(dopriParams, 'MinStep'), defaultOpts.MinStep = dopriParams.MinStep; end
    if isfield(integratorOptions, 'verbose'), verbose = integratorOptions.verbose; end

    % Override defaults with user-provided integratorOptions
    currentOpts = defaultOpts;
    if nargin > 4 && ~isempty(integratorOptions)
        user_fields = fieldnames(integratorOptions);
        for i = 1:length(user_fields)
            currentOpts.(user_fields{i}) = integratorOptions.(user_fields{i});
        end
    end
    
    if verbose
        fprintf('  Adaptive DoPri8(7) Integration Setup (Returning Steps):\n');
        fprintf('    Range: [%.3f, %.3f]\n', z_start, z_end);
        fprintf('    Particle Initial KE (E0): %.4f. Using U = q*V_elec - E0.\n', E0);
        fprintf('    DA Order: %d, nVars (state): %d\n', daOrder, nVars);
        fprintf('    RelTol: %.2e, AbsTol: %.2e\n', currentOpts.RelTol, currentOpts.AbsTol);
    end
    % Print other options if they exist

    % --- Initialization ---
    r_da = DiffAlg.var(1, daOrder, nVars); % r (DA var 1)
    pr_da = DiffAlg.var(2, daOrder, nVars);% p_r = dr/dz (DA var 2)
    S_initial = {r_da, pr_da}; % Initial state (identity map)

    % Package shared parameters for the derivative function
    derivParams.q_charge = q_charge;
    derivParams.E0 = E0;
    derivParams.daOrder = daOrder;
    derivParams.nVars = nVars;
    derivParams.bemParams = bemParams;
    derivParams.bemTable = bemTable;

    % Define the derivative function handle
    f_ode = @(z_eval, S_eval) state_deriv_dopri(z_eval, S_eval, derivParams);

    % --- Perform Full Adaptive Integration ---
    tic;
    if verbose
        fprintf('    Starting adaptive integration (will return steps)...\n');
    end

    % Call the modified integrator that returns steps
    [~, z_steps, S_steps, stats] = integrate_dopri87_da(f_ode, z_span, S_initial, currentOpts);
    % S_final is implicitly S_steps{end}
    % z_final is implicitly z_steps(end)

    integration_time = toc;
    if verbose
        fprintf('    Finished adaptive integration in %.2f seconds.\n', integration_time);
        fprintf('    Integration ended at z = %.6f (Target z_end = %.6f)\n', z_steps(end), z_end);
        fprintf('    Returned %d intermediate steps.\n', numel(z_steps));
        if abs(z_steps(end) - z_end) > abs(z_end)*1e-9
             warning('Integrator did not reach the specified z_end accurately.');
        end
    end

    % --- Check for NaNs in the final state (optional but good) ---
    S_final = S_steps{end};
    if any(isnan(S_final{1}.Series)) || any(isnan(S_final{2}.Series))
        error('NaN detected in the final state S_final (S_steps{end}). Integration failed.');
    end

    % No need to extract R_total here, as the primary outputs are the steps.
    % If R_total is needed, it can be extracted from S_steps{end} outside this function.
    
    if verbose
        fprintf('  Returning z_steps and S_steps for ray tracing.\n');
    end

end % --- End of main function ---


% ======================================================
% --- Local Helper Function for State Derivative ---
% ======================================================
% IMPORTANT: Keep this function definition identical to the one in the
% previous modification of run_dopri87_da_integration.m or ensure it's
% accessible via MATLAB's path.

function dS = state_deriv_dopri(z_eval, S_eval, params)
    % ... (Keep the exact same implementation as provided before) ...
    % Calculates the state derivative {dr/dz, dp/dz} at a specific z_eval.
    % Unpack parameters
    q_charge = params.q_charge;
    E0 = params.E0;
    daOrder = params.daOrder;
    nVars = params.nVars; % Should be 2
    qVs = params.bemParams.qVs;
    r_center_bem = params.bemTable.r_center;
    z_center_bem = params.bemTable.z_center;

    % Unpack state variables (DA objects, nVars=2)
    % r = S_eval{1}; % Not needed directly if evaluating DA objects below
    p = S_eval{2}; % p = dr/dz

    % 1. Calculate Axial Coefficients for V_elec at z_eval
    Nmax_coeffs = max([daOrder * 2, 10, 1]);
    Coeffs_eval = axial_potential_coeffs(Nmax_coeffs, z_eval, qVs, r_center_bem, z_center_bem);

    % 2. Construct DA V_elec(r) at z_eval
    r_da_local = DiffAlg.var(1, daOrder, nVars);
    V_elec_at_z_da = DiffAlg.one(nVars, daOrder) * Coeffs_eval(1);
    r_sq_local = r_da_local * r_da_local;
    r_pow_2k_local = r_sq_local;
    max_k_pot = floor(daOrder / 2);
    for k = 1:max_k_pot
        if (2*k + 1) > numel(Coeffs_eval), break; end
        if isnan(Coeffs_eval(2*k + 1)), continue; end
        V_2k = factorial(2*k) * Coeffs_eval(2*k + 1);
        coeff_k = ((-1)^k) / ( (4^k) * (factorial(k)^2) );
        if ~isfinite(coeff_k) || ~isfinite(V_2k), continue; end
        term_da = r_pow_2k_local * (coeff_k * V_2k);
        V_elec_at_z_da = V_elec_at_z_da + term_da;
        if k < max_k_pot && (2*(k+1) <= daOrder)
            r_pow_2k_local = r_pow_2k_local * r_sq_local;
        end
    end

    % 3. Construct DA U(r) at z_eval
    U_at_z_da = (V_elec_at_z_da * q_charge) - E0 * DiffAlg.one(nVars, daOrder);

    % 4. Calculate DA dU/dr at z_eval
    Ur_at_z_da = U_at_z_da.differentiate(1);

    % 5. Calculate SCALAR dU/dz at z_eval
    if numel(Coeffs_eval) >= 2 && ~isnan(Coeffs_eval(2))
        V_elec_z_eval = Coeffs_eval(2);
    else
        V_elec_z_eval = 0;
    end
    Uz_at_z_scalar = q_charge * V_elec_z_eval;

    % --- Now compute the derivatives using the current state S_eval ---
    one = DiffAlg.one(nVars, daOrder);

    drdz = p; % Eq 1

    one_plus_p2 = one + p * p;
    U_eval = U_at_z_da;
    Ur_eval = Ur_at_z_da;

    U0_const_eval = U_eval.Series(1);
    nan_da = []; % Initialize flag/value for NaN return
    if U0_const_eval >= -1e-12 || abs(U_eval.Series(1)) < 1e-15
        % Potential issue: reflection or division by zero
        if U0_const_eval >= -1e-12
             % warning('state_deriv_dopri:potential', 'U >= 0 (%.2e) at z=%.3f. Possible reflection.', U0_const_eval, z_eval);
        end
         if abs(U_eval.Series(1)) < 1e-15
             % warning('state_deriv_dopri:denominator', 'Denominator (2*U) near zero (%.2e) at z=%.3f. Division unstable.', U_eval.Series(1), z_eval);
         end
         nan_series = NaN(size(drdz.Series));
         nan_da = DiffAlg(nan_series, daOrder, nVars);
         dS = {nan_da, nan_da};
         return;
    end

    term_in_parentheses = Ur_eval - (p * Uz_at_z_scalar);
    numerator_da = one_plus_p2 * term_in_parentheses;
    denominator_da = U_eval * 2.0;
    dpdz = numerator_da / denominator_da; % Eq 2

    % Check if dpdz contains NaN/Inf after division
     if any(~isfinite(dpdz.Series))
         warning('state_deriv_dopri:nonfinite_dpdz', 'Non-finite value detected in dpdz at z=%.4f.', z_eval);
         nan_series = NaN(size(drdz.Series));
         nan_da = DiffAlg(nan_series, daOrder, nVars);
         dS = {nan_da, nan_da}; % Return NaN state derivative
         return;
     end

    dS = {drdz, dpdz};

end % End state_deriv_dopri