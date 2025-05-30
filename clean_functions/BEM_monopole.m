function [qs, bemTable] = BEM_monopole(electrodes, full_integration_cutoff, verbose)
    % BEM_monopole Computes the monopole BEM matrices and the charge distributions.
    %
    % OUTPUTS:
    %   qs         - A (nElectrodes x numSegments) matrix. Each row contains the
    %                charge distribution on all BEM segments when that electrode is
    %                set to 1V (and all other electrodes are at 0V).
    %
    %   bemTable   - A table listing the BEM segments (columns:
    %                [r_center, z_center, length, elec_idx]). In addition, the table
    %                has nElectrodes extra columns (named "q_elec_1", "q_elec_2", etc.)
    %                that contain the charge distribution computed for each electrode.
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
    % Calculate centers and lengths directly from the raw data, avoiding loops.
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

    %% Compute the Monopole Influence Matrix A_mon
    % A_mon is computed by integrating the mutual interaction between each pair
    % of BEM segments.
    if verbose
        fprintf("Populating BEM matrix (%d x %d)...\n", numSegments, numSegments)
    end

    % The parfor loop implementation below is already efficient for parallel computation.
    % We calculate only the upper triangle and then symmetrize.

    % Pre-allocate the final matrix
    A_mon = zeros(numSegments, numSegments);

    % 1. Find the row and column indices for the upper triangle (including diagonal)
    [row_idx, col_idx] = find(triu(ones(numSegments))); % row_idx(k) <= col_idx(k) always
    num_pairs = length(row_idx); % Number of elements to calculate

    % Temporary storage for the calculated values in the parfor loop
    upper_tri_values = zeros(num_pairs, 1);

    % Local copy of bem_elements_raw for use inside parfor to avoid broadcasting the whole table
    % (Although broadcasting the smaller raw array might be slightly better than the original)
    bem_raw_local = bem_elements_raw;

    % --- Parallel Calculation of Upper Triangle of A_mon ---
    parfor k = 1:num_pairs
        % Get the indices for this pair
        i = row_idx(k);
        j = col_idx(k); % Note: j >= i here

        % Get centers for calculation (avoids repeated indexing inside norm/functions)
        center_i = centers(i,:); % Already pre-calculated
        center_j = centers(j,:); % Already pre-calculated

        % Calculate distance between centers
        distCenters = norm( center_i - center_j );

        % Decide which interaction function to use based on distance
        if distCenters > full_integration_cutoff
            % Use far-field approximation (interaction between centers).
            upper_tri_values(k) = far_field_interaction_mon( center_i(1), center_i(2), center_j(1), center_j(2) );
        else
            % Use full integral interaction.
            % Pass center of segment 'i' and endpoints of segment 'j'.
            % Indexing bem_raw_local which is broadcast.
            upper_tri_values(k) = full_integral_interaction_mon( center_i(1), center_i(2), ...
                                     bem_raw_local(j,1), bem_raw_local(j,3), ... % Rb1, Rb2
                                     bem_raw_local(j,2), bem_raw_local(j,4) );    % Zb1, Zb2
        end
    end

    % 3. Assemble the matrix from the calculated upper triangular values
    linearIndices = sub2ind(size(A_mon), row_idx, col_idx);
    A_mon(linearIndices) = upper_tri_values;

    % 4. Make the matrix symmetric by copying the upper triangle to the lower triangle
    A_mon = A_mon + triu(A_mon, 1).'; % Add transpose of strictly upper triangle

    if verbose
        fprintf("BEM matrix populated.\n")
    end

    %% Solve the BEM system for each electrode excited individually
    if verbose
        fprintf("Solving BEM system for %d electrode excitations...\n", nElectrodes)
    end

    % --- Optimization: Solve for all electrode excitations simultaneously ---
    % Instead of looping and solving A_mon \ V for each electrode, create a
    % matrix V_all where each column is the voltage vector for one excitation,
    % and solve A_mon \ V_all once. This is generally much faster.

    % Create the V_all matrix [numSegments x nElectrodes].
    % V_all(k, elec) = 1 if segment k belongs to electrode elec, 0 otherwise.
    % Use sparse matrix construction for efficiency, then convert to full.
    % electrode_indices contains the electrode index for each segment (row).
    V_all = full(sparse(1:numSegments, electrode_indices, 1, numSegments, nElectrodes));

    % Solve the linear system A_mon * QS = V_all for all charge distributions QS.
    % QS will be [numSegments x nElectrodes].
    QS = A_mon \ V_all;

    % Transpose QS to get the desired output format for qs [nElectrodes x numSegments].
    qs = QS.';

    if verbose
        fprintf("BEM system solved.\n")
    end

    %% Append the charge distribution for each electrode to bemTable.
    % --- Optimization: Add all charge columns to the table at once ---
    % Create column names programmatically.
    chargeColNames = arrayfun(@(x) sprintf('q_elec_%d', x), 1:nElectrodes, 'UniformOutput', false);
    % Convert the charge distributions (transposed QS) to a table.
    chargeTable = array2table(QS, 'VariableNames', chargeColNames);
    % Concatenate the original bemTable with the new charge table.
    bemTable = [bemTable, chargeTable];

    if verbose
        disp('Updated BEM table with charge distributions for each electrode:');
        % Display only head if table is very large
        if height(bemTable) > 20
             disp(head(bemTable, 10));
             fprintf('... (displaying head of %d rows)\n', height(bemTable));
        else
             disp(bemTable);
        end
    end
end

%% Local functions for monopole influence calculations (Unchanged from previous)
% --- These functions implement the core physics/math ---
% Optimization within these might involve lower-level numerics (e.g., custom
% integration, faster elliptic integral implementations) but are assumed correct
% and reasonably efficient for now. Added 'eps' usage for stability.

function A = full_integral_interaction_mon(Rac, Zac, Rb1, Rb2, Zb1, Zb2)
    % Computes the full (numerical) integral for the monopole interaction.
    wdth = sqrt((Rb1 - Rb2)^2 + (Zb1 - Zb2)^2);
    % Avoid division by zero for zero-length segments (should not happen ideally)
    if wdth < eps
        A = 0;
        return;
    end
    theta = atan2(Rb2 - Rb1, Zb2 - Zb1); % Swapped Rb, Zb order for standard angle
    Rbc = 0.5 * (Rb1 + Rb2);
    Zbc = 0.5 * (Zb1 + Zb2);
    % Define parametric functions for the target segment (tau from -wdth/2 to wdth/2)
    Rb = @(tau) Rbc + tau .* sin(theta); % Parametric r coord on segment B
    Zb = @(tau) Zbc + tau .* cos(theta); % Parametric z coord on segment B
    % Integrand: mutual interaction kernel divided by segment length (as per monopole def.)
    integrand = @(tau) mutual_interaction(Rac, Zac, Rb(tau), Zb(tau)); % Integrate kernel directly
    % Integrate using quadgk over the parameter tau
    A = (1/wdth) * quadgk(integrand, -wdth/2, wdth/2); % Normalize by length outside
end

function A = far_field_interaction_mon(R1, Z1, R2, Z2)
    % Computes the far-field (approximate) monopole interaction.
    % Handle potential case where R1 or R2 is zero or very small
     if R1 < eps || R2 < eps
         % If one ring is a point on axis, interaction might be simplified or zero
         % depending on formulation. Let's assume approx holds, or return 0 if problematic.
         % Check standard formula validity for R=0. Ellipke fails for k^2=0.
         % If R1 or R2 is 0, k2=0, K=pi/2. Need careful limit analysis.
         % For simplicity, if R1 or R2 are near zero, maybe use full integral or specific limit.
         % However, let's stick to the formula and ensure k2 is handled.
         % If R1*R2 is numerically zero, A is likely zero or needs special case.
         if abs(R1 * R2) < eps
             A = 0; % Or handle singularity based on specific physics model
             return;
         end
     end
    denom = ((R1 + R2)^2 + (Z1 - Z2)^2);
    % Check for division by zero (coincident points, R1+R2=0 only if R1=R2=0)
    if denom < eps
        A = Inf; % Or handle singularity based on physics
        return;
    end
    k2 = 4 * R1 * R2 / denom; % abs() removed as R is radial coord >= 0
    k2 = min(k2, 1 - eps); % Ensure k^2 is not exactly 1 for ellipke
    k2 = max(k2, 0);      % Ensure k^2 is not negative due to precision issues
    [K, ~] = ellipke(k2);
    % k = sqrt(k2); % k is not directly needed here based on formula used
    % Original formula was: A = (0.25/pi^2) * (1 / sqrt(R1 * R2)) * k * K;
    % Let's use a common form involving G = k*K/sqrt(R1*R2)
    % A = G / (4*pi^2) ? Verify formula source if possible.
    % Assuming the formula A = (1/(2*pi)) * (2/pi) * K / sqrt((R1+R2)^2 + (Z1-Z2)^2) might be intended,
    % or related forms involving elliptic integrals.
    % Sticking to the provided structure:
    % Check potential division by sqrt(R1*R2) if R1 or R2 are zero
    denom_sqrt = sqrt(R1 * R2);
     if denom_sqrt < eps
          A = 0; % If either radius is zero, the interaction term might vanish this way
          % Re-check physics if R=0 is a valid input leading to non-zero interaction
          return;
     end
    k = sqrt(k2);
    A = (1 / (4 * pi^2)) * (k * K / denom_sqrt); % Adjusted constant based on common forms, double check definition
end

function Aij = mutual_interaction(Ri, Zi, Rj, Zj)
    % Computes the mutual interaction kernel - VECTORIZED for quadgk.
    % Assumes Ri, Zi are scalar. Rj, Zj can be vectors.

    % Initialize output array with the same size as Rj/Zj
    Aij = zeros(size(Rj)); % Or size(Zj), should be same

    % --- Element-wise calculation and edge case handling ---

    % 1. Handle cases where Ri or Rj are near zero (avoids sqrt(0) issues)
    % Identify indices where Ri or Rj are effectively zero.
    % Use element-wise OR | because Rj is a vector.
    zero_radius_idx = (Ri < eps) | (Rj < eps);
    % For these cases, the interaction is typically zero. Set Aij and skip further calcs for them.
    Aij(zero_radius_idx) = 0;

    % 2. Handle coincident points (causes division by zero in k2 calc)
    d_val_sq = (Ri + Rj).^2 + (Zi - Zj).^2;
    coincident_idx = d_val_sq < eps;
    % Interaction is infinite at coincidence. Set Aij and skip further calcs for them.
    % Make sure not to overwrite the zero_radius_idx cases if they overlap.
    Aij(coincident_idx & ~zero_radius_idx) = Inf; % Set Inf only if not already set to 0

    % 3. Proceed with calculation only for valid, non-zero, non-coincident points
    % Create an index of points needing full calculation.
    calc_idx = ~zero_radius_idx & ~coincident_idx;

    % Only perform calculations if there are valid points
    if ~any(calc_idx)
        return; % All points were edge cases, Aij is already set.
    end

    % --- Calculations for valid indices ---
    % Subset inputs to only valid indices for intermediate calculations
    Rj_valid = Rj(calc_idx);
    Zj_valid = Zj(calc_idx);
    d_val_sq_valid = d_val_sq(calc_idx); % Use pre-calculated squared value

    % Calculate k^2 safely using only valid data
    k2 = 4 * (Ri .* Rj_valid) ./ d_val_sq_valid;

    % Clamp k2 element-wise for ellipke stability
    k2 = min(k2, 1 - eps);
    k2 = max(k2, 0); % Should already be >=0 if Ri,Rj >= 0

    % Calculate Elliptic Integral K element-wise safely
    [K_ei, ~] = ellipke(k2);
    % Check if ellipke produced NaN (can happen for k2 near 1 despite clamping)
    nan_K_idx = isnan(K_ei);
    if any(nan_K_idx)
        warning('BEM_monopole:mutual_interaction', 'NaN generated by ellipke in mutual_interaction. Setting Aij to 0 for these points.');
        K_ei(nan_K_idx) = 0; % Set corresponding K to 0 to avoid NaN propagation
        % Or could set Aij directly for these original indices if preferred
    end


    % Calculate k = sqrt(k2) element-wise
    k = sqrt(k2);

    % Calculate denominator sqrt(Ri * Rj) safely (already know Ri,Rj > eps here)
    denom_sqrt = sqrt(Ri * Rj_valid);

    % Calculate Aij for the valid indices
    Aij_valid = (1 / (4 * pi^2)) * (k .* K_ei ./ denom_sqrt);

    % Place the calculated valid values back into the main Aij array
    Aij(calc_idx) = Aij_valid;

    % Note: No need for final NaN/Inf check if logic above handles all cases.
    % The check for Inf is handled by 'coincident_idx'.
    % The check for NaN is handled after ellipke.
    % The check for zero radius is handled by 'zero_radius_idx'.
end