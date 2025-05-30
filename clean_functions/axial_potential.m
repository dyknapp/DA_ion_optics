function U = axial_potential(zs, qVs, bem_rs, bem_zs)
    % Calculates axial potential U at points zs due to BEM charges qVs.
    % Uses a loop over zs points for memory efficiency.
    arguments
        zs (1, :) double       % Axial positions (row vector, 1xM)
        qVs (:, 1) double      % Charge on each BEM segment (column vector, Nx1)
        bem_rs (:, 1) double   % Radius of each BEM segment (column vector, Nx1)
        bem_zs (:, 1) double   % Z-position of each BEM segment (column vector, Nx1)
    end

    num_bem = size(bem_zs, 1); % Number of BEM segments (N)
    num_zs = numel(zs);        % Number of axial points (M)

    % --- Input Validation ---
    if size(bem_rs, 1) ~= num_bem || size(qVs, 1) ~= num_bem
        error('Mismatch in number of BEM elements (rows) for qVs, bem_rs, bem_zs.');
    end
    if isempty(zs) || isempty(qVs)
        U = zeros(1, num_zs); % Return empty or appropriately sized zero vector
        return;
    end

    % --- Memory-Efficient Calculation ---
    % Pre-allocate output vector
    U = zeros(1, num_zs);

    % Pre-calculate square of BEM radii (constant in the loop)
    bem_rs_sq = bem_rs.^2; % Nx1 vector

    % Loop over each axial point where potential is requested
    for j = 1:num_zs
        % Current axial point
        z_point = zs(j);

        % Calculate (z_point - bem_zs)^2 for all segments element-wise
        delta_z_sq = (z_point - bem_zs).^2; % Nx1 vector

        % Calculate denominator = pi * sqrt( (z_point-bem_zs)^2 + bem_rs^2 )
        denom = 4 * pi * sqrt(delta_z_sq + bem_rs_sq); % Nx1 vector

        % Avoid potential division by zero (e.g., z_point coincides with a
        % zero-radius ring center). Replace near-zero denominators with a small number.
        denom(denom < eps) = eps;

        % Calculate potential contributions from all rings to this z_point
        contributions = qVs ./ denom; % Nx1 vector

        % Sum contributions from all N segments to get total potential at zs(j)
        U(j) = sum(contributions);
    end
    % U is now a 1xM row vector containing the potential at each zs point.
end