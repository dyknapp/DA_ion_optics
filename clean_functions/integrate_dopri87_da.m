function [S_final, z_steps, S_steps, stats] = integrate_dopri87_da(f_deriv, z_span, S_initial, options)
% Integrates a system of ODEs for DA objects using Dormand-Prince 8(7)
% with adaptive step size control. Stores intermediate steps.
% INPUTS:
%   f_deriv   : Derivative function handle f(z, S) -> {dS1/dz, dS2/dz,...}
%   z_span    : [z_start, z_end] integration interval
%   S_initial : Initial state cell array {S1_da, S2_da, ...}
%   options   : Struct with integrator options (RelTol, AbsTol, InitialStep, MaxStep, MinStep, Stats)
% OUTPUTS:
%   S_final   : Final state cell array at z_end
%   z_steps   : Vector [z_start, z_step1, z_step2, ..., z_final] of z values at each accepted step.
%   S_steps   : Cell array {S_initial, S_step1, S_step2, ..., S_final} of DA state maps at each z_step.
%   stats     : Struct with integration statistics (if requested)

    % --- Default Options ---
    opts = struct();
    opts.RelTol = 1e-6;
    opts.AbsTol = 1e-9;
    opts.MaxStep = abs(z_span(2) - z_span(1)) / 10;
    opts.MinStep = abs(z_span(2) - z_span(1)) * 1e-12;
    opts.InitialStep = [];
    opts.Stats = 'off';
    opts.MaxRetries = 5;
    opts.SafetyFactor = 0.9;
    opts.MinScaleFactor = 0.2;
    opts.MaxScaleFactor = 5.0;
    % Option for pre-allocating steps (can improve performance slightly)
    opts.PreallocateSteps = 1000; % Estimate number of steps

    % Override defaults with user options
    if nargin > 3 && ~isempty(options)
        user_fields = fieldnames(options);
        for i = 1:length(user_fields)
            opts.(user_fields{i}) = options.(user_fields{i});
        end
    end

    % --- Initialization ---
    z_start = z_span(1);
    z_end = z_span(2);
    S = S_initial; % Current state map
    z = z_start;   % Current z value
    h_dir = sign(z_end - z_start);
    h_max = abs(opts.MaxStep) * h_dir;
    h_min = abs(opts.MinStep) * h_dir;

    % --- Initialize Statistics and Output Storage ---
    stats.nSteps = 0;
    stats.nFailed = 0;
    stats.nFunEvals = 0;

    % Pre-allocate or initialize storage for intermediate steps
    if opts.PreallocateSteps > 0
        z_steps = zeros(opts.PreallocateSteps, 1);
        S_steps = cell(opts.PreallocateSteps, 1);
        storage_capacity = opts.PreallocateSteps;
    else
        z_steps = []; % Grow dynamically
        S_steps = {}; % Grow dynamically
        storage_capacity = 0;
    end
    step_count = 1; % Start storing from index 1
    z_steps(step_count) = z_start;
    S_steps{step_count} = S_initial;


    % Get coefficients once
    [c, A, b, b_star] = get_dopri87_coeffs(); % Assuming this function is available
    coeffs.c = c; coeffs.A = A; coeffs.b = b; coeffs.b_star = b_star;
    nStages = 13; % Number of stages in DoPri8(7)

    % --- Automatic Initial Step Size ---
    if isempty(opts.InitialStep)
        try
            k1_deriv = f_deriv(z, S); stats.nFunEvals = stats.nFunEvals + 1;
             if any(isnan(k1_deriv{1}.Series)) || any(isnan(k1_deriv{2}.Series))
                 error('NaN detected during initial derivative evaluation.');
             end
            scale = opts.AbsTol + max_abs_coeffs(S) * opts.RelTol;
            if scale == 0, scale = opts.AbsTol; end
            err_deriv_norm = max_abs_coeffs(k1_deriv);
            if err_deriv_norm < 1e-5, d1 = 1e-6; else, d1 = err_deriv_norm; end
            h0 = min(abs(h_max), abs(z_end-z)*0.1);
            h = min(abs(h_max), h0 * (opts.AbsTol / d1)^(1/8)); % Rough guess
            h = h * h_dir;
            if strcmp(opts.Stats, 'on'), fprintf('Automatic initial step size: %.4e\n', h); end
        catch ME_init_step
             warning('Automatic initial step size calculation failed: %s. Using MaxStep/10.', ME_init_step.message);
             h = h_max / 10.0;
        end
    else
        h = opts.InitialStep * h_dir;
    end
    % Ensure initial step is within bounds
    h = max(abs(h), abs(h_min)) * h_dir;
    h = min(abs(h), abs(h_max)) * h_dir;
    % Prevent zero initial step if z_start = z_end
    if abs(h) < eps && abs(z_start - z_end) > eps
        h = h_min * h_dir * 10; % Use a small step if calc resulted in zero
        h = min(abs(h), abs(h_max)) * h_dir;
    end

    % --- Integration Loop ---
    integration_complete = false;
    while ~integration_complete
        % Check if we are very close to the end
        if abs(z - z_end) < abs(h_min) * 0.5
            integration_complete = true;
            break; % Reached end point within tolerance
        end

        % Ensure step doesn't overshoot z_end significantly
        if (z + h - z_end) * h_dir > 0 % If current step will overshoot
            h = z_end - z; % Adjust step to land exactly on z_end
            if abs(h) < eps % If remaining distance is effectively zero
                 integration_complete = true;
                 break;
            end
        end
        % Prevent step size from becoming too small unless near the end
        if abs(h) < abs(h_min)
            if abs(z - z_end) > abs(h_min) * 1.1
                warning('Integrator step size h=%.2e reached minimum allowed h_min=%.2e at z=%.4f before reaching z_end=%.4f.', abs(h), abs(h_min), z, z_end);
                h = h_min * h_dir; % Force minimum step
                if (z + h - z_end) * h_dir > 0, h = z_end - z; end % Adjust if min step overshoots
                if abs(h) < eps, integration_complete = true; break; end
            else
                h = z_end - z; % Take final small step to reach end
                if abs(h) < eps, integration_complete = true; break; end
            end
        end

        % --- Attempt Step ---
        step_accepted = false;
        retries = 0;
        S_next = S; % Initialize S_next in case of loop failure

        while ~step_accepted && retries < opts.MaxRetries
            % Perform one DoPri8(7) step calculation
            [S_try, S_err] = dopri87_step_da(z, h, S, f_deriv, coeffs);
            stats.nFunEvals = stats.nFunEvals + nStages;

            % Check for NaNs from step function (indicates f_deriv failed)
             if any(isnan(S_try{1}.Series)) || any(isnan(S_try{2}.Series))
                 warning('NaN detected in integrator step evaluation at z=%.4f. Reducing step size.', z);
                 h = h * opts.MinScaleFactor * 0.5; % Drastic reduction
                 h = max(abs(h), abs(h_min)) * h_dir; % Keep within bounds
                 stats.nFailed = stats.nFailed + 1;
                 retries = retries + 1;
                 if retries >= opts.MaxRetries
                     error('Integrator failed: NaN encountered and max retries reached at z=%.4f.', z);
                 end
                 continue; % Retry with smaller h
             end

            % Calculate scalar error norm
            error_norm = calculate_da_error_norm(S_err, S_try, opts.AbsTol, opts.RelTol);

            % Estimate optimal step size for NEXT step
            if error_norm < eps % Avoid division by zero / large scaling if error is tiny
                h_new_abs = abs(h) * opts.MaxScaleFactor;
            else
                scale_factor = opts.SafetyFactor * (1.0 / error_norm)^(1/8);
                % Limit scaling factor
                scale_factor = max(opts.MinScaleFactor, min(opts.MaxScaleFactor, scale_factor));
                h_new_abs = abs(h) * scale_factor;
            end
            h_new_abs = max(abs(h_min), min(abs(h_max), h_new_abs)); % Clamp between min/max
            h_new = h_new_abs * h_dir;

            % Check if step is accepted
            if error_norm <= 1.0 % Step successful
                step_accepted = true;
                stats.nSteps = stats.nSteps + 1;
                z = z + h;      % Update z
                S = S_try;      % Update state
                h = h_new;      % Use optimal size for next step

                % Store the accepted step results
                step_count = step_count + 1;
                % Check if storage needs resizing (if not preallocated or exceeded)
                if step_count > storage_capacity
                    if opts.PreallocateSteps > 0
                        % Double the storage size
                        new_capacity = storage_capacity * 2;
                         fprintf('Resizing storage from %d to %d steps at z=%.4f\n', storage_capacity, new_capacity, z);
                        z_steps_new = zeros(new_capacity, 1);
                        S_steps_new = cell(new_capacity, 1);
                        z_steps_new(1:storage_capacity) = z_steps;
                        S_steps_new(1:storage_capacity) = S_steps;
                        z_steps = z_steps_new;
                        S_steps = S_steps_new;
                        storage_capacity = new_capacity;
                    else
                        % Allow dynamic growth (can be slow)
                    end
                end
                z_steps(step_count) = z;
                S_steps{step_count} = S; % Store the cumulative map at this z

                if strcmp(opts.Stats, 'on') && mod(stats.nSteps, 50) == 0 % Print progress periodically
                    fprintf(' Step %d: z=%.6f, h=%.4e, err=%.3e\n', stats.nSteps, z, h, error_norm);
                end

            else % Step failed, reduce step size and retry
                stats.nFailed = stats.nFailed + 1;
                retries = retries + 1;
                h_old = h;
                h = h_new; % Reduce step size for retry (already scaled down)
                if strcmp(opts.Stats, 'on')
                   fprintf(' Step failed at z=%.6f (h=%.2e). err=%.3e > 1.0. Retrying with h=%.4e (Retry %d/%d)\n', z, h_old, error_norm, h, retries, opts.MaxRetries);
                end
                if abs(h) <= abs(h_min) * 1.01 && error_norm > 1.0
                    warning('Cannot achieve requested tolerance with minimum step size at z=%.4f. Error=%.3e.', z, error_norm);
                    if retries >= opts.MaxRetries
                        error('Integrator failed: Max retries reached at minimum step size.');
                    end
                end
            end % End if error_norm <= 1.0
        end % End while ~step_accepted

        if ~step_accepted % If loop finished due to max retries
            error('Integrator failed: Step not accepted after %d retries at z=%.4f.', opts.MaxRetries, z);
        end

    end % End while ~integration_complete

    % --- Finalize Outputs ---
    S_final = S; % The last accepted state

    % Trim unused preallocated storage
    if step_count < storage_capacity
        z_steps = z_steps(1:step_count);
        S_steps = S_steps(1:step_count);
    end

    % Ensure the final point is exactly z_end if we stopped slightly short
    if abs(z_steps(end) - z_end) > eps && abs(z_steps(end) - z_end) < abs(h_min) * 0.6
        % If the loop terminated just before z_end because the remainder was too small,
        % take one final exact step IF the state derivative is still valid at z_end.
        h_final = z_end - z_steps(end);
        try
            [S_final_exact, ~] = dopri87_step_da(z_steps(end), h_final, S_steps{end}, f_deriv, coeffs);
             if ~any(isnan(S_final_exact{1}.Series)) && ~any(isnan(S_final_exact{2}.Series))
                S_final = S_final_exact;
                z_steps(end) = z_end; % Correct final z
                S_steps{end} = S_final; % Update final state map
                if strcmp(opts.Stats, 'on'), fprintf('Took final small step h=%.2e to reach z_end exactly.\n', h_final); end
             else
                 warning('NaN detected trying final step to z_end. Returning state at z=%.6f', z_steps(end));
             end
        catch ME_final_step
             warning('Error during final step to z_end: %s. Returning state at z=%.6f', ME_final_step.message, z_steps(end));
        end
    elseif abs(z_steps(end)-z_end) > eps
         warning('Final z (%.6f) differs significantly from z_end (%.6f)', z_steps(end), z_end);
    end


    if strcmp(opts.Stats, 'on')
        fprintf('Integration finished.\n');
        fprintf('  Steps: %d successful, %d failed.\n', stats.nSteps, stats.nFailed);
        fprintf('  Function Evaluations: %d\n', stats.nFunEvals);
        fprintf('  Final z: %.8f\n', z_steps(end));
    end

end

% --- Helper: Calculate Max Abs Coeff for DA objects/cell array ---
function maxval = max_abs_coeffs(S_da)
    if iscell(S_da)
        maxval = 0;
        for i = 1:length(S_da)
            % Check if Series is empty or DA object is invalid before accessing
            if ~isempty(S_da{i}) && isprop(S_da{i}, 'Series') && ~isempty(S_da{i}.Series)
                 maxval = max(maxval, max(abs(S_da{i}.Series)));
            end
        end
    elseif isa(S_da, 'DiffAlg') && ~isempty(S_da.Series)
        maxval = max(abs(S_da.Series));
    else
        maxval = 0; % Or throw error
    end
     if ~isfinite(maxval), maxval = 0; end % Handle potential Inf/NaN in Series
end

% --- Helper: Calculate Error Norm ---
function err_norm = calculate_da_error_norm(S_err, S_next, abs_tol, rel_tol)
    max_err_scaled = 0;
    for i = 1:length(S_err) % Loop through state variables (r, pr)
        % Use max_abs_coeffs which handles empty/invalid cases
        scale_i = abs_tol + rel_tol * max_abs_coeffs(S_next{i});
        if scale_i < eps, scale_i = abs_tol; end % Avoid division by zero or near-zero scale

        % Check S_err{i} validity
        if ~isempty(S_err{i}) && isprop(S_err{i}, 'Series') && ~isempty(S_err{i}.Series)
            err_coeffs_i = abs(S_err{i}.Series);
            % Check for non-finite error coefficients
             if any(~isfinite(err_coeffs_i))
                 warning('Non-finite values found in error estimate S_err{%d}. Returning Inf norm.', i);
                 max_err_scaled = Inf;
                 break; % Exit loop if Inf error encountered
             end
            err_scaled_coeffs_i = err_coeffs_i / scale_i;
            max_err_scaled = max(max_err_scaled, max(err_scaled_coeffs_i));
        end
    end
    err_norm = max_err_scaled;
     if ~isfinite(err_norm), err_norm = Inf; end % Ensure Inf is returned if loop broke early
end

% --- State Addition/Scaling Helpers (Needed by dopri87_step_da) ---
% Include state_add and state_scale functions here if dopri87_step_da doesn't define them itself.
% function Ssum = state_add(S1, S2) ... end
% function Sscaled = state_scale(S, factor) ... end

% Need state_add and state_scale functions (copy from previous script)
function Ssum = state_add(S1, S2)
    Ssum = cell(size(S1));
    for i = 1:numel(S1)
        Ssum{i} = S1{i} + S2{i};
    end
end
function Sscaled = state_scale(S, factor)
     Sscaled = cell(size(S));
    for i = 1:numel(S)
        Sscaled{i} = S{i} * factor;
    end
end