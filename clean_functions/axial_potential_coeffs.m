function C_n = axial_potential_coeffs(Nmax, z0, qVs, bem_rs, bem_zs)
% Calculates coefficients C_n for the axial potential series expansion
% U(z) = sum_{n=0}^{Nmax} C_n * (z - z0)^n
% based on BEM segment charges and positions.

arguments
    Nmax (1,1) {mustBeInteger, mustBeNonnegative} % Max order n of the expansion
    z0 (1,1) double          % Expansion point on the z-axis
    qVs (:, 1) double         % Charge on each BEM segment (Nx1 column vector)
    bem_rs (:, 1) double      % Radius of each BEM segment (Nx1 column vector)
    bem_zs (:, 1) double      % Z-position of each BEM segment (Nx1 column vector)
end

num_bem = size(bem_zs, 1); % Number of BEM segments (N)

% --- Input Validation ---
if size(bem_rs, 1) ~= num_bem || size(qVs, 1) ~= num_bem
    error('Mismatch in number of BEM elements (rows) for qVs, bem_rs, bem_zs.');
end
if Nmax < 0
    error('Nmax must be non-negative.');
end

% Initialize total coefficients C_n (for n=0 to Nmax) -> size Nmax+1
C_n_total = zeros(Nmax + 1, 1);

% --- Loop through each BEM segment ---
for i = 1:num_bem
    R_i = bem_rs(i);
    bemZ_i = bem_zs(i);
    qVs_i = qVs(i);

    % Calculate Y[n] sequence for this segment up to Nmax
    Yn_i = calculate_Yn_sequence(R_i, bemZ_i, z0, Nmax);

    % Calculate c_n = Y[n] / (4*pi) for this segment
    % Apply the 1/(4*pi) scaling factor
    cn_i = Yn_i / (4 * pi);

    % Add the charge-weighted coefficients of this segment to the total
    C_n_total = C_n_total + qVs_i * cn_i;
end

C_n = C_n_total; % Return coefficients as a column vector (index k corresponds to n=k-1)

end

% --- Helper function to compute Y[n] sequence ---
function Yn = calculate_Yn_sequence(R, bemZ, z0, Nmax)
    % Calculates Y[n] for n=0 to Nmax using the recurrence relation derived
    % from the Mathematica DifferenceRoot object.

    % Pre-allocate Yn vector (MATLAB index k corresponds to n = k-1)
    Yn = zeros(Nmax + 1, 1);

    % Calculate intermediate variables
    dz = bemZ - z0;
    d2 = dz^2 + R^2;

    % --- Handle potential singularity ---
    if abs(d2) < 1e-30 % Avoid division by zero if expansion point z0
                      % coincides with the center of a zero-radius ring (R=0).
        warning('Potential singularity detected (d2=dz^2+R^2 is near zero) for R=%.3g, z=%.3g, z0=%.3g. Returning NaN.', R, bemZ, z0);
        Yn(:) = NaN;
        return;
    end

    sqrt_d2 = sqrt(d2);

    % --- Initial Conditions (Y[0], Y[1]) ---
    % Y[0] corresponds to Yn(1)
    Yn(1) = 1 / sqrt_d2;

    % Y[1] corresponds to Yn(2) - only calculate if Nmax >= 1
    if Nmax >= 1
        Yn(2) = dz / (d2 * sqrt_d2); % or dz / (d2^(3/2))
    end

    % --- Recurrence Relation Loop ---
    % Loop to calculate Y[n+2] (index n+3) for n = 0, 1, ..., Nmax-2
    for n = 0:Nmax-2
        % Y[n]   is Yn(n+1)
        % Y[n+1] is Yn(n+2)
        % Y[n+2] is Yn(n+3)

        numerator = (3 + 2*n) * dz * Yn(n+2) - (1 + n) * Yn(n+1);

        % Denominator is (2 + n) * d2. Check for safety although d2 != 0.
        denominator = (2 + n) * d2;
        if abs(denominator) < 1e-30
             warning('Recurrence relation denominator is near zero for n=%d. R=%.3g, z=%.3g, z0=%.3g. Returning NaN.', n, R, bemZ, z0);
             Yn(n+3:end) = NaN; % Set remaining terms to NaN
             return; % Stop calculation for this ring
        end

        Yn(n+3) = numerator / denominator;
    end
end