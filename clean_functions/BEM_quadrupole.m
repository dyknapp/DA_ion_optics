% BEM_quadrupole.m

function [qs_quad, bemTable] = BEM_quadrupole(electrodes, full_integration_cutoff, verbose)
    % BEM_quadrupole Computes the quadrupole BEM matrices and the quadrupole source distributions.
    % Based on the structure of BEM_monopole.m
    %
    % OUTPUTS:
    %   qs_quad    - A (nElectrodes x numSegments) matrix. Each row contains the
    %                quadrupole source distribution on all BEM segments when that
    %                electrode is set to 1V (and all other electrodes are at 0V).
    %
    %   bemTable   - A table listing the BEM segments (columns:
    %                [r_center, z_center, length, elec_idx]). In addition, the table
    %                has nElectrodes extra columns (named "q_quad_elec_1", etc.)
    %                that contain the quadrupole source distribution computed for each electrode.
    %
    % INPUTS:
    %   electrodes - A 1xN vector of Electrode objects. Each electrode object is expected to
    %                have properties "rs" (a vector of radial positions) and "segments"
    %                (an N_seg-by-4 matrix whose columns are: [r_start, z_start, r_end, z_end]).
    %   full_integration_cutoff - A scalar distance for switching between full integration
    %                and far-field approximations (default = 1.0).
    %   verbose    - (Optional) Logical flag to display intermediate results. Default is true.
    %
    arguments
        electrodes (1,:) Electrode
        full_integration_cutoff (1, 1) double = 1.0
        verbose logical = true
    end

    nElectrodes = numel(electrodes);

    %% Build the BEM elements array and calculate segment properties
    % Collect all segments from all electrodes and calculate geometric properties.
    % We pre-calculate total segments for efficient allocation.
    totalSegments = 0;
    for idx = 1:nElectrodes
        % Assuming electrode.segments is [N x 4] where N = numel(rs)-1
        totalSegments = totalSegments + size(electrodes(idx).segments, 1);
    end

    % Pre-allocate arrays for segment data
    bem_elements_raw = zeros(totalSegments, 4); % [r_start, z_start, r_end, z_end]
    electrode_indices = zeros(totalSegments, 1); % [elec_idx]
    rowCounter = 0;

    % Populate the raw segment data and electrode indices
    for idx = 1:nElectrodes
        segs = electrodes(idx).segments; % Expected size: [numSegs x 4]
        numSegs = size(segs, 1);
        rows = (rowCounter+1):(rowCounter+numSegs);
        bem_elements_raw(rows, :) = segs;
        electrode_indices(rows) = idx;  % Mark these segments as belonging to electrode idx
        rowCounter = rowCounter + numSegs;
    end

    % --- Optimization: Vectorized calculation of segment centers and lengths ---
    centers = [0.5 * (bem_elements_raw(:, 1) + bem_elements_raw(:, 3)), ... % r_center
               0.5 * (bem_elements_raw(:, 2) + bem_elements_raw(:, 4))];   % z_center
    lengths = sqrt((bem_elements_raw(:, 1) - bem_elements_raw(:, 3)).^2 + ...
                   (bem_elements_raw(:, 2) - bem_elements_raw(:, 4)).^2);

    numSegments = totalSegments; % Use the pre-calculated total

    % --- Modification: Create bemTable with center, length, and index ---
    bemTable = table(centers(:, 1), centers(:, 2), lengths, electrode_indices, ...
                     'VariableNames', {'r_center', 'z_center', 'length', 'elec_idx'});

    if verbose
        disp('BEM elements (r_center, z_center, length, elec_idx):');
        disp(bemTable);
    end

    %% Compute the Quadrupole Influence Matrix A_quad
    % A_quad is computed by integrating the mutual quadrupole interaction between
    % each pair of BEM segments.
    if verbose
        fprintf("Populating Quadrupole BEM matrix (%d x %d)...\n", numSegments, numSegments)
    end

    % Pre-allocate the final matrix
    A_quad = zeros(numSegments, numSegments);

    % 1. Find the row and column indices for the upper triangle (including diagonal)
    [row_idx, col_idx] = find(triu(ones(numSegments))); % row_idx(k) <= col_idx(k) always
    num_pairs = length(row_idx); % Number of elements to calculate

    % Temporary storage for the calculated values in the parfor loop
    upper_tri_values = zeros(num_pairs, 1);

    % Local copy of bem_elements_raw for use inside parfor
    bem_raw_local = bem_elements_raw;

    % --- Parallel Calculation of Upper Triangle of A_quad ---
    parfor k = 1:num_pairs
        % Get the indices for this pair
        i = row_idx(k);
        j = col_idx(k); % Note: j >= i here

        % Get centers for calculation
        center_i = centers(i,:);
        center_j = centers(j,:);

        % Calculate distance between centers
        distCenters = norm( center_i - center_j );

        % Decide which interaction function to use based on distance
        if distCenters > full_integration_cutoff
            % Use far-field approximation (interaction between centers).
            % The far_field_interaction function returns [A_mon, A_quad].
            [~, A_q_far] = far_field_interaction( center_i(1), center_i(2), center_j(1), center_j(2) );
            upper_tri_values(k) = A_q_far;
        else
            % Use full integral interaction for quadrupole.
            % Pass center of segment 'i' and endpoints of segment 'j'.
            upper_tri_values(k) = full_integral_interaction_quad( center_i(1), center_i(2), ...
                                     bem_raw_local(j,1), bem_raw_local(j,3), ... % Rb1, Rb2
                                     bem_raw_local(j,2), bem_raw_local(j,4) );    % Zb1, Zb2
        end
    end

    % 3. Assemble the matrix from the calculated upper triangular values
    linearIndices = sub2ind(size(A_quad), row_idx, col_idx);
    A_quad(linearIndices) = upper_tri_values;

    % 4. Make the matrix symmetric by copying the upper triangle to the lower triangle
    A_quad = A_quad + triu(A_quad, 1).'; % Add transpose of strictly upper triangle

    if verbose
        fprintf("Quadrupole BEM matrix populated.\n")
    end

    %% Solve the BEM system for each electrode excited individually
    if verbose
        fprintf("Solving Quadrupole BEM system for %d electrode excitations...\n", nElectrodes)
    end

    % Create the V_all matrix [numSegments x nElectrodes].
    V_all = full(sparse(1:numSegments, electrode_indices, 1, numSegments, nElectrodes));

    % Solve the linear system A_quad * QS = V_all for all quadrupole source distributions QS.
    % QS will be [numSegments x nElectrodes].
    QS_quad = A_quad \ V_all;

    % Transpose QS_quad to get the desired output format [nElectrodes x numSegments].
    qs_quad = QS_quad.';

    if verbose
        fprintf("Quadrupole BEM system solved.\n")
    end

    %% Append the quadrupole source distribution for each electrode to bemTable.
    % Create column names programmatically.
    quadSourceColNames = arrayfun(@(x) sprintf('q_quad_elec_%d', x), 1:nElectrodes, 'UniformOutput', false);
    % Convert the quadrupole source distributions (transposed QS_quad) to a table.
    quadSourceTable = array2table(QS_quad, 'VariableNames', quadSourceColNames);
    % Concatenate the original bemTable with the new quadrupole source table.
    bemTable = [bemTable, quadSourceTable];

    if verbose
        disp('Updated BEM table with quadrupole source distributions for each electrode:');
        % Display only head if table is very large
        if height(bemTable) > 20
             disp(head(bemTable, 10));
             fprintf('... (displaying head of %d rows)\n', height(bemTable));
        else
             disp(bemTable);
        end
    end
end

%% Local functions for quadrupole influence calculations (Provided by user)
% --- These functions implement the core physics/math for quadrupole ---

function A = full_integral_interaction_quad(Rac, Zac, Rb1, Rb2, Zb1, Zb2)
    % Full integration for quadrupole interaction when segments are close.
    % Computes the integral of the quadrupole ring interaction kernel over segment B.
    wdth = sqrt((Rb1 - Rb2)^2 + (Zb1 - Zb2)^2);
    % Avoid division by zero for zero-length segments
    if wdth < eps
        A = 0;
        return;
    end
    theta = atan2(Rb2 - Rb1, Zb2 - Zb1); % Angle of segment B
    Rbc = 0.5*(Rb1 + Rb2); % Center r of segment B
    Zbc = 0.5*(Zb1 + Zb2); % Center z of segment B

    % Parameterize segment B using tau from -wdth/2 to wdth/2
    Rb = @(tau) Rbc + tau.*sin(theta); % Use tau directly, not 0.5*tau
    Zb = @(tau) Zbc + tau.*cos(theta); % Use tau directly, not 0.5*tau

    % Define the integrand: quadrupole ring interaction kernel
    integrand = @(tau) quadrupole_ring_interaction(Rac, Zac, Rb(tau), Zb(tau));

    % Integrate using quadgk over the parameter tau and normalize by length
    try
        integral_value = quadgk(integrand, -wdth/2, wdth/2);
        A = (1/wdth) * integral_value; % Normalize by length AFTER integration
    catch ME
        warning('BEM_quadrupole:quadgkError', 'quadgk failed in full_integral_interaction_quad: %s. Returning 0.', ME.message);
        A = 0;
    end
end

function [A_mon, A_quad] = far_field_interaction(R1, Z1, R2, Z2)
    % Computes both the monopole and quadrupole interactions for far-field segments.
    % Based on interaction between segment centers (R1, Z1) and (R2, Z2).

    % Handle potential case where R1 or R2 is zero or very small
     if R1 < eps || R2 < eps
         if abs(R1 * R2) < eps
             A_mon = 0;
             A_quad = 0;
             return;
         end
     end

    denom = ( (R1 + R2).^2 + (Z1 - Z2).^2 );
     if denom < eps % Check for coincident points
        warning('BEM_quadrupole:far_field_interaction:Coincidence', 'Far-field denominator near zero. R1=%.2e, R2=%.2e, Z1=%.2e, Z2=%.2e. Returning NaN.', R1, R2, Z1, Z2);
        A_mon = NaN;
        A_quad = NaN;
        return;
    end

    k2 = 4 * abs(R1 .* R2) ./ denom; % abs() used in original snippet
    k2 = min(k2, 1 - eps); % Ensure k^2 is not exactly 1 for ellipke
    k2 = max(k2, 0);      % Ensure k^2 is not negative

    % Handle case k2 -> 0 safely for Quadrupole part prefactor
    k2_quad_denom = k2.^(3/2);
    if k2 < eps
        % If k2 is zero (e.g., R1 or R2 is zero), avoid division by zero.
        % The quadrupole interaction should go to zero in this limit.
        A_quad = 0;
    else
        % Quadrupole part calculation (only if k2 > eps)
        try
            [eK, eE] = ellipke(k2);
             % Check for NaN/Inf from ellipke
             if any(isnan(eK)) || any(isnan(eE)) || any(isinf(eK)) || any(isinf(eE))
                 warning('BEM_quadrupole:far_field_interaction:EllipkeNaN', 'ellipke returned NaN/Inf for k2=%.2e. R1=%.2e, R2=%.2e. Returning 0 for A_quad.', k2, R1, R2);
                 A_quad = 0;
             else
                prefactor_quad = pi ./ (3 * sqrt(abs(R1 .* R2)) .* k2_quad_denom + eps); % Add eps for sqrt safety
                term1_quad = 16 * (-2 + k2) .* eE;
                term2_quad = 2  * (16 - 16*k2 + 3*k2.^2) .* eK;
                A_quad = prefactor_quad .* (term1_quad + term2_quad);
                 % Final check for NaN/Inf in A_quad
                if isnan(A_quad) || isinf(A_quad)
                     warning('BEM_quadrupole:far_field_interaction:ResultNaN', 'A_quad calculation resulted in NaN/Inf. R1=%.2e, R2=%.2e. Returning 0.', R1, R2);
                     A_quad = 0;
                end
             end
        catch ME_ellip
            warning('BEM_quadrupole:far_field_interaction:EllipkeError', 'Error in ellipke: %s. R1=%.2e, R2=%.2e. Returning 0 for A_quad.', ME_ellip.message, R1, R2);
            A_quad = 0;
        end
    end

    % Monopole part calculation (can use the one from BEM_monopole if preferred)
    % Re-calculating here based on the provided function structure.
     denom_sqrt_mon = sqrt(abs(R1 .* R2));
     if denom_sqrt_mon < eps
         A_mon = 0;
     else
        try
            if k2 > eps % Recalculate eK if not already done or if k2 was 0
                 [eK, ~] = ellipke(k2);
            elseif ~exist('eK','var') % Case k2=0
                 eK = pi/2;
            end
             k = sqrt(k2);
             A_mon = (1 / (4 * pi^2)) * (k * eK / denom_sqrt_mon); % Matched to BEM_monopole version
              % Final check for NaN/Inf in A_mon
             if isnan(A_mon) || isinf(A_mon)
                 warning('BEM_quadrupole:far_field_interaction:ResultNaNMon', 'A_mon calculation resulted in NaN/Inf. R1=%.2e, R2=%.2e. Returning 0.', R1, R2);
                 A_mon = 0;
             end
        catch ME_ellip_mon
            warning('BEM_quadrupole:far_field_interaction:EllipkeErrorMon', 'Error in monopole ellipke: %s. R1=%.2e, R2=%.2e. Returning 0 for A_mon.', ME_ellip_mon.message, R1, R2);
            A_mon = 0;
        end
     end
end

function A = quadrupole_ring_interaction(R1, Z1, R2, Z2)
    % Quadrupole potential factor kernel between two rings.
    % Assumes R1, Z1 scalar; R2, Z2 can be vectors from quadgk.
    A = zeros(size(R2)); % Initialize output

    % --- Element-wise calculation for vector R2/Z2 ---
    zero_radius_idx = (R1 < eps) | (R2 < eps);
    A(zero_radius_idx) = 0;

    denom_sq = ((R1 + R2).^2 + (Z1 - Z2).^2);
    coincident_idx = denom_sq < eps;
    % Set A=NaN for coincident points (or Inf, but NaN might be safer for quadgk)
    A(coincident_idx & ~zero_radius_idx) = NaN;

    calc_idx = ~zero_radius_idx & ~coincident_idx;
    if ~any(calc_idx), return; end % Exit if all points are edge cases

    % Subset inputs for calculation
    R2_valid = R2(calc_idx);
    Z2_valid = Z2(calc_idx);
    denom_sq_valid = denom_sq(calc_idx);

    % Calculate k^2 safely
    k2 = abs(4 * (R1 .* R2_valid) ./ denom_sq_valid); % Use abs() as in original
    k2 = min(k2, 1 - eps);
    k2 = max(k2, 0);

    A_valid = zeros(size(R2_valid)); % Initialize temporary result

    % --- Calculate for non-zero k2 ---
    k2_nonzero_idx = k2 > eps;
    if any(k2_nonzero_idx)
        k2_nz = k2(k2_nonzero_idx);
        R2_nz = R2_valid(k2_nonzero_idx); % Need R2 corresponding to non-zero k2
        try
            [Kval, Eval] = ellipke(k2_nz);
             if any(isnan(Kval)) || any(isnan(Eval)) || any(isinf(Kval)) || any(isinf(Eval))
                warning('BEM_quadrupole:quad_ring:EllipkeNaN', 'ellipke returned NaN/Inf for k2=%.2e. R1=%.2e, R2=%.2e. Setting A to 0 for these.', k2_nz, R1, R2_nz);
                % Find original indices corresponding to NaN/Inf and set A directly? More complex.
                % For now, set corresponding A_valid elements to 0.
                A_valid_nz = zeros(size(Kval)); % Default to 0
                nan_inf_ellip_idx = isnan(Kval) | isnan(Eval) | isinf(Kval) | isinf(Eval);
             else
                k2_quad_denom_nz = k2_nz.^(3/2);
                denom_sqrt_nz = sqrt(abs(R1 .* R2_nz)); % Use R2_nz
                prefactor = pi ./ (3 * denom_sqrt_nz .* k2_quad_denom_nz + eps); % Add eps for safety
                term1 = 16 * (-2 + k2_nz) .* Eval;
                term2 = 2  * (16 - 16*k2_nz + 3*k2_nz.^2) .* Kval;
                A_valid_nz = prefactor .* (term1 + term2);
                % Check for NaN/Inf result
                 nan_inf_A_idx = isnan(A_valid_nz) | isinf(A_valid_nz);
                 if any(nan_inf_A_idx)
                    warning('BEM_quadrupole:quad_ring:ResultNaN', 'A calculation resulted in NaN/Inf. R1=%.2e, R2=%.2e. Setting to 0.', R1, R2_nz(nan_inf_A_idx));
                    A_valid_nz(nan_inf_A_idx) = 0;
                 end
             end
             A_valid(k2_nonzero_idx) = A_valid_nz; % Place results back

        catch ME_ellip_ring
             warning('BEM_quadrupole:quad_ring:EllipkeError', 'Error in ellipke: %s. R1=%.2e, R2=%.2e. Setting A to 0.', ME_ellip_ring.message, R1, R2_nz);
             A_valid(k2_nonzero_idx) = 0; % Set corresponding elements to 0
        end
    end
    % A_valid already contains 0 for k2==0 cases (implicitly).

    % Place valid results back into the main A array
    A(calc_idx) = A_valid;

    % --- Unstable for very small R: just shunt it to zero ---
    % (Handled by zero_radius_idx check at the beginning)
    % A(abs(R1.*R2) < 1e-6) = 0; % Redundant if zero_radius_idx is used.

end