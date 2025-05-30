function [r_positions_combined, z_positions_combined, initial_states_out] = ...
    trace_ray_fan_multisegment_da_map(traceParams, integrationResults1, ...
    S_maps_bender_segmented, z_coords_bender, integrationResults2, ...
    visParams, max_eval_order, show)
%trace_ray_fan_multisegment_da_map Traces a fan of rays through three segments:
%   1. DA integration result 1
%   2. A sequence of DA maps (e.g., for a bender)
%   3. DA integration result 2
%   Plots rays color-coded by rank (Sort: Descending r0, then Descending s0).
%
%   INPUTS:
%       traceParams: Struct with ray tracing parameters:
%                    .num_r_rays, .num_s_rays, .max_r0_dev, .max_s0_dev,
%                    .center_r0, .center_s0 (p_r = dr/dz)
%       integrationResults1: Struct with integrator output for first segment:
%                            .z_steps (vector), .S_steps (cell array of maps {r,p})
%       S_maps_bender_segmented: Cell array of DA map pairs {{r_map1,p_map1}, {r_map2,p_map2},...}
%                                for the bender segments.
%       z_coords_bender: Vector of z-coordinates at the boundaries of the bender segments.
%                        Should have length = numel(S_maps_bender_segmented) + 1.
%       integrationResults2: Struct with integrator output for third segment:
%                            .z_steps (vector), .S_steps (cell array of maps {r,p})
%       visParams: Struct with visualization info (electrodes_full_stack, simParams)
%       max_eval_order: (Optional) Maximum DA order for map evaluation.
%       show: (Optional) Boolean, true to display figure (default=true)
%
%   OUTPUTS:
%       r_positions_combined: Matrix (num_z_points_total x num_rays_total) of radial positions.
%       z_positions_combined: Vector (num_z_points_total x 1) of corresponding z coordinates.
%       initial_states_out:   Matrix (2 x num_rays_total) of the initial [r0; p0] used.
%
%   Displays Figure 2 with the ray fan plot.
%
%   Dependencies: DiffAlg class (with evaluate method).

    % --- Extract Parameters ---
    num_r_rays = traceParams.num_r_rays;
    num_s_rays = traceParams.num_s_rays;
    max_r0_dev = traceParams.max_r0_dev;
    max_s0_dev = traceParams.max_s0_dev;
    center_r0 = traceParams.center_r0;
    center_s0 = traceParams.center_s0; % Note: s0 is p_r = dr/dz
    z_steps1 = integrationResults1.z_steps(:); % Ensure column
    S_steps1 = integrationResults1.S_steps; % Cell array {r_map, pr_map} cumulative from start
    num_bender_segments = numel(S_maps_bender_segmented);
    z_coords_bender = z_coords_bender(:); % Ensure column
    z_steps2 = integrationResults2.z_steps(:); % Ensure column
    S_steps2 = integrationResults2.S_steps; % Cell array {r_map, pr_map} cumulative from z_bend_end
    % Use the full electrode stack for visualization boundaries
    electrodes = visParams.electrodes_full_stack;
    simParams = visParams.simParams;

    % --- Input Validation & Defaults ---
    map_order = S_steps1{1}{1}.Order; % Get DA order from first map
    if nargin < 8
        show = true;
    end
    if nargin < 7 || isempty(max_eval_order)
         max_eval_order = map_order; % Default to full map order
         if show, fprintf('  Using default evaluation order: %d (full map order)\n', max_eval_order); end
    elseif max_eval_order > map_order
         warning('trace_ray_fan_multisegment:OrderExceeded', ...
                 'Requested evaluation order (%d) exceeds map order (%d). Using map order.', max_eval_order, map_order);
         max_eval_order = map_order;
    elseif max_eval_order < 0
          warning('trace_ray_fan_multisegment:NegativeOrder', ...
                 'Requested evaluation order (%d) is negative. Using order 0.', max_eval_order);
         max_eval_order = 0;
     else
         if show, fprintf('  Using specified evaluation order: %d\n', max_eval_order); end
     end

     % Add checks for DA order consistency across all segments if needed
     if numel(z_coords_bender) ~= (num_bender_segments + 1)
         error('trace_ray_fan_multisegment:BenderCoordMismatch', ...
               'Number of z_coords_bender (%d) must be n_segments (%d) + 1.', numel(z_coords_bender), num_bender_segments);
     end
     z_bend_start = z_coords_bender(1);
     z_bend_end = z_coords_bender(end);
     % Optional check for boundary alignment
     % if abs(z_steps1(end) - z_bend_start) > 1e-6 || abs(z_steps2(1) - z_bend_end) > 1e-6
     %      warning('trace_ray_fan_multisegment:BoundaryMismatch', ...
     %              'Z-boundary mismatch between integration results and bender start/end points (z1_end=%.6f, z_bend_start=%.6f; z2_start=%.6f, z_bend_end=%.6f).', ...
     %              z_steps1(end), z_bend_start, z_steps2(1), z_bend_end);
     % end

    % --- Define Initial Ray Fan Grid ---
    if num_r_rays <= 1, r_devs = 0; else r_devs = linspace(-max_r0_dev, max_r0_dev, num_r_rays); end
    if num_s_rays <= 1, s_devs = 0; else s_devs = linspace(-max_s0_dev, max_s0_dev, num_s_rays); end
    [R_dev_grid, S_dev_grid] = meshgrid(r_devs, s_devs);
    initial_r0s = center_r0 + R_dev_grid(:);
    initial_s0s = center_s0 + S_dev_grid(:); % s0 = p_r0
    initial_states = [initial_r0s'; initial_s0s']; % 2 x num_rays_total
    initial_states_out = initial_states;
    num_rays_total = size(initial_states, 2);
    if show, fprintf('  Generating initial fan: %d x %d = %d total rays.\n', num_r_rays, num_s_rays, num_rays_total); end
    if show, fprintf('    r0 range: [%.3f, %.3f] (Center: %.3f)\n', min(initial_r0s), max(initial_r0s), center_r0); end
    if show, fprintf('    p_r0 range: [%.3f, %.3f] (Center: %.3f)\n', min(initial_s0s), max(initial_s0s), center_s0); end

    % --- Trace Rays by Evaluating DA Maps Sequentially ---
    num_z_points1 = numel(z_steps1);
    num_z_points_bender = num_bender_segments + 1; % Number of boundaries
    num_z_points2 = numel(z_steps2); % Includes start point z_bend_end

    z_positions_combined = [z_steps1; z_coords_bender(2:end); z_steps2(2:end)];
    num_z_points_total = numel(z_positions_combined);
    r_positions_combined = zeros(num_z_points_total, num_rays_total);
    if show, fprintf('  Tracing %d rays by evaluating DA maps (Order <= %d) through %d total points...\n', num_rays_total, max_eval_order, num_z_points_total); end
    tic;
    calculation_errors = false;
    X_current_state = zeros(2, num_rays_total);

    % --- Trace Section 1 ---
    if show, fprintf('    Tracing section 1 (%d points)...\n', num_z_points1); end
    r_positions_combined(1,:) = initial_states(1,:);
    for k = 1:num_rays_total
        initial_state_k = initial_states(:, k)';
        for i = 2:num_z_points1
             if any(isnan(r_positions_combined(i:end, k))), continue; end
             S_map = S_steps1{i};
             try
                 r_new = S_map{1}.evaluate(initial_state_k, max_eval_order);
                 p_new = S_map{2}.evaluate(initial_state_k, max_eval_order);
                 if ~isfinite(r_new) || ~isfinite(p_new)
                     warning('trace_ray_fan_multisegment:EvalNaNInf1', 'NaN/Inf result for ray %d at z=%.4f (Sec 1). Setting rest to NaN.', k, z_steps1(i));
                     r_positions_combined(i:end, k) = NaN;
                     calculation_errors = true; continue;
                 end
                 r_positions_combined(i, k) = r_new;
                 if i == num_z_points1
                     X_current_state(:, k) = [r_new; p_new];
                 end
             catch ME_eval
                 fprintf('Error eval DA map (Sec 1) ray %d at z=%.4f: %s\n', k, z_steps1(i), ME_eval.message);
                 r_positions_combined(i:end, k) = NaN;
                 calculation_errors = true;
             end
        end % i loop
    end % k loop

    % --- Trace Bender Section ---
    if show, fprintf('    Tracing bender section (%d segments)...\n', num_bender_segments); end
    start_idx_bender_in_combined = num_z_points1;
    for i = 1:num_bender_segments
        global_idx_end_of_segment = start_idx_bender_in_combined + i;
        S_map_i = S_maps_bender_segmented{i};
        for k = 1:num_rays_total
             if any(isnan(r_positions_combined(global_idx_end_of_segment:end, k))), continue; end
             state_at_segment_start_k = X_current_state(:, k)';
             try
                 r_new = S_map_i{1}.evaluate(state_at_segment_start_k, max_eval_order);
                 p_new = S_map_i{2}.evaluate(state_at_segment_start_k, max_eval_order);
                 if ~isfinite(r_new) || ~isfinite(p_new)
                     warning('trace_ray_fan_multisegment:EvalNaNInfBender', 'NaN/Inf result for ray %d at z=%.4f (Bender seg %d). Setting rest to NaN.', k, z_coords_bender(i+1), i);
                     r_positions_combined(global_idx_end_of_segment:end, k) = NaN;
                     calculation_errors = true; continue;
                 end
                 r_positions_combined(global_idx_end_of_segment, k) = r_new;
                 X_current_state(:, k) = [r_new; p_new];
             catch ME_eval
                 fprintf('Error eval DA map (Bender seg %d) ray %d at z=%.4f: %s\n', i, k, z_coords_bender(i+1), ME_eval.message);
                 r_positions_combined(global_idx_end_of_segment:end, k) = NaN;
                 calculation_errors = true;
             end
        end % k loop
    end % i loop

    % --- Trace Section 3 ---
    if show, fprintf('    Tracing section 3 (%d points)...\n', num_z_points2 - 1); end
    start_idx_sec3_in_combined = start_idx_bender_in_combined + num_bender_segments;
    state_at_sec3_start = X_current_state;
    for k = 1:num_rays_total
        state_after_bend_k = state_at_sec3_start(:, k)';
        for i = 2:num_z_points2
            global_idx = start_idx_sec3_in_combined + i - 1;
             if any(isnan(r_positions_combined(global_idx:end, k))), continue; end
            S_map = S_steps2{i};
            try
                r_new = S_map{1}.evaluate(state_after_bend_k, max_eval_order);
                if ~isfinite(r_new)
                     warning('trace_ray_fan_multisegment:EvalNaNInf2', 'NaN/Inf result for ray %d at z=%.4f (Sec 3). Setting rest to NaN.', k, z_steps2(i));
                     r_positions_combined(global_idx:end, k) = NaN;
                     calculation_errors = true; continue;
                 end
                r_positions_combined(global_idx, k) = r_new;
            catch ME_eval
                 fprintf('Error eval DA map (Sec 3) ray %d at z=%.4f: %s\n', k, z_steps2(i), ME_eval.message);
                 r_positions_combined(global_idx:end, k) = NaN;
                 calculation_errors = true;
            end
        end % i loop
    end % k loop
    trace_time = toc;
    if show, fprintf('  Finished ray tracing in %.2f seconds.\n', trace_time); end
     if calculation_errors, fprintf('  NOTE: Errors occurred during DA map evaluation (NaN/Inf results).\n'); end

    % --- Plot Ray Fan ---
    if show
        % figure(2); clf;
        hold on;

        % === MODIFIED SECTION START: Rank-based Coloring ===
        % Define sorting criteria: Descending r0, then descending s0
        % Use negative values for descending sort with sortrows
        sort_criteria = [-initial_states(1,:)', -initial_states(2,:)']; % Columns: [-r0, -s0]
        original_indices = (1:num_rays_total)';

        % Sort based on criteria and get the order of original indices
        % sorted_map gives the original indices in the desired sorted order
        [~, sorted_map] = sortrows(sort_criteria, [1 2]);

        % Assign colors based on the rank (position in the sorted order)
        cmap = cool(256); % Or choose another colormap (e.g., parula, viridis)
        num_colors = size(cmap, 1);
        color_indices = zeros(1, num_rays_total); % Preallocate vector for color indices

        for rank = 1:num_rays_total
            original_index = sorted_map(rank); % Get the original index of the ray at this rank

            % Normalize the rank to [0, 1] for colormap indexing
            if num_rays_total <= 1
                 normalized_rank = 0.5;
            else
                 normalized_rank = (rank - 1) / (num_rays_total - 1); % Rank 1 -> 0, Rank N -> 1
            end

            % Determine the color index for this ray based on its rank
            color_indices(original_index) = round(normalized_rank * (num_colors - 1)) + 1;
        end
        % Now color_indices(k) holds the colormap index for the k-th ray according to its rank
        % === MODIFIED SECTION END ===

        % Plot rays individually using the rank-based color index
        for k = 1:num_rays_total
            ray_r_positions = r_positions_combined(:, k);
            valid_indices = ~isnan(ray_r_positions);
            if any(valid_indices)
                plot(z_positions_combined(valid_indices), ray_r_positions(valid_indices), ...
                     '-', 'Color', cmap(color_indices(k), :), 'LineWidth', 0.5);
            end
        end

        % Plot electrode geometry
        if isfield(visParams, 'electrodes_full_stack') && ~isempty(visParams.electrodes_full_stack)
            electrodes_vis = visParams.electrodes_full_stack;
        else
            electrodes_vis = visParams.electrodes;
            warning('Using electrode geometry from visParams.electrodes. Ensure this includes bender representation.');
        end
        for k = 1:numel(electrodes_vis)
            elec = electrodes_vis(k);
            if isprop(elec,'rs') && isprop(elec,'zs') && ~isempty(elec.rs) && ~isempty(elec.zs)
                plot(elec.zs, elec.rs, '-k', 'LineWidth', 1.5);
                plot(elec.zs, -elec.rs, '-k', 'LineWidth', 1.5);
            end
        end

        xlabel('Axial Position z [mm]');
        ylabel('Radial Position r [mm]');
        title('Ray Fan For Ion Trap Injection');
        grid on;
        axis tight;

        ylim_curr = ylim;
        ylim_max_abs = max(abs(ylim_curr));
        ylim_required_extent = ylim_max_abs;
        if isfield(simParams, 'R') && ~isempty(simParams.R)
             ylim_required_extent = max(ylim_required_extent, simParams.R * 1.1);
        end
         ylim([-ylim_required_extent, ylim_required_extent]);

        xline(z_positions_combined(1), ':k', {'Start'}, 'LabelVerticalAlignment', 'bottom');
        xline(z_bend_start, '--r', {'Bender', 'Start'}, 'LabelVerticalAlignment', 'bottom');
        xline(z_bend_end, '--r', {'Bender', 'End'}, 'LabelVerticalAlignment', 'bottom');
        xline(z_positions_combined(end), ':k', {'End'}, 'LabelVerticalAlignment', 'bottom');

        hold off;
        fprintf('  Multi-segment ray fan plot generated in Figure 2 (colored by rank: high r0, then high s0).\n'); % Updated description
    end
end