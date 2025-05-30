function [S_next, S_err] = dopri87_step_da(z, h, S, f_deriv, coeffs)
% Performs one step of the Dormand-Prince 8(7) method for DA objects.
% INPUTS:
%   z        : Current independent variable value (z_start_step)
%   h        : Step size
%   S        : Current state cell array {r_da, pr_da}
%   f_deriv  : Derivative function handle f(z, S) -> {dr/dz, dp/dz}
%   coeffs   : Struct containing DoPri8(7) coefficients (c, A, b, b_star)
% OUTPUTS:
%   S_next   : State cell array at z+h using 8th order result
%   S_err    : Estimated error cell array (8th order - 7th order result)

    c = coeffs.c;
    A = coeffs.A;
    b = coeffs.b;
    b_star = coeffs.b_star;
    s = 13; % Number of stages

    % Store intermediate derivative results (k_i) scaled by h
    K = cell(s, 1); % Each K{i} will be {h*dr/dz_i, h*dp/dz_i}

    % Stage 1
    try
        k_deriv_1 = f_deriv(z, S);
    catch ME
        fprintf('Error in f_deriv (k1) at z=%.4f: %s\n', z, ME.message); rethrow(ME);
    end
    if any(isnan(k_deriv_1{1}.Series)) || any(isnan(k_deriv_1{2}.Series))
        S_next = {k_deriv_1{1},k_deriv_1{2}}; S_err = S_next; return; % Propagate NaN
    end
    K{1} = state_scale(k_deriv_1, h);

    % Stages 2 to s=13
    for i = 2:s
        % Calculate state S_temp = S + sum(A(i,j)*K{j}) for j=1 to i-1
        S_temp = S; % Start with initial state for this stage
        for j = 1:(i-1)
            if A(i,j) ~= 0
                S_temp = state_add(S_temp, state_scale(K{j}, A(i,j)));
            end
        end

        % Evaluate derivative at stage i
        z_eval = z + c(i) * h;
        try
            k_deriv_i = f_deriv(z_eval, S_temp);
        catch ME
             fprintf('Error in f_deriv (k%d) at z=%.4f: %s\n', i, z_eval, ME.message); rethrow(ME);
        end
        if any(isnan(k_deriv_i{1}.Series)) || any(isnan(k_deriv_i{2}.Series))
             S_next = {k_deriv_i{1},k_deriv_i{2}}; S_err = S_next; return; % Propagate NaN
        end
        K{i} = state_scale(k_deriv_i, h);
    end

    % Calculate final 8th order result: S_next = S + sum(b(i)*K{i})
    S_next = S;
    for i = 1:s
        if b(i) ~= 0
            S_next = state_add(S_next, state_scale(K{i}, b(i)));
        end
    end

    % Calculate 7th order result for error estimation: S_7 = S + sum(b_star(i)*K{i})
    S_7 = S;
     for i = 1:s
        if b_star(i) ~= 0
            S_7 = state_add(S_7, state_scale(K{i}, b_star(i)));
        end
    end

    % Error estimate: S_err = S_next - S_7
    S_err = state_add(S_next, state_scale(S_7, -1.0));

end

% --- State Addition Helper ---
function Ssum = state_add(S1, S2)
    Ssum = {S1{1} + S2{1}, S1{2} + S2{2}};
end

% --- State Scaling Helper ---
function Sscaled = state_scale(S, factor)
    Sscaled = {S{1} * factor, S{2} * factor};
end