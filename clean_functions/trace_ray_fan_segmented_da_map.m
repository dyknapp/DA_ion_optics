function [r_positions_combined, z_positions_combined, initial_states_out] = trace_ray_fan_segmented_da_map(traceParams, integrationResults1, S_map_bender, integrationResults2, z_bend_start, z_bend_end, visParams, max_eval_order)
%trace_ray_fan_segmented_da_map Traces a fan of rays using segmented DA maps evaluated sequentially.
%   Plots rays color-coded by initial radial position.
%
%   INPUTS:
%       traceParams: Struct with ray tracing parameters:
%                    .num_r_rays, .num_s_rays, .max_r0_dev, .max_s0_dev,
%                    .center_r0, .center_s0 (p_r = dr/dz)
%       integrationResults1: Struct with integrator output for first segment:
%                            .z_steps (vector), .S_steps (cell array of maps)
%       S_map_bender: Cell array {r_bender_da, p_bender_da} DA map for bender section.
%       integrationResults2: Struct with integrator output for third segment:
%                            .z_steps (vector), .S_steps (cell array of maps)
%       z_bend_start: Z-coordinate at the start of the bender.
%       z_bend_end:   Z-coordinate at the end of the bender.
%       visParams: Struct with visualization info (electrodes_full_stack, simParams)
%       max_eval_order: (Optional) Maximum DA order for map evaluation.
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

    z_steps1 = integrationResults1.z_steps;
    S_steps1 = integrationResults1.S_steps; % Cell array {r_map, pr_map} cumulative from start
    z_steps2 = integrationResults2.z_steps;
    S_steps2 = integrationResults2.S_steps; % Cell array {r_map, pr_map} cumulative from z_bend_end

    % Use the full electrode stack for visualization boundaries
    electrodes = visParams.electrodes_full_stack;
    simParams = visParams.simParams;

    % --- Input Validation & Defaults ---
    map_order = S_steps1{1}{1}.Order; % Get DA order from maps
    if nargin < 8 || isempty(max_eval_order)
         max_eval_order = map_order; % Default to full map order
         fprintf('  Using default evaluation order: %d (full map order)\n', max_eval_order);
    elseif max_eval_order > map_order
         warning('trace_ray_fan_segmented_da_map:OrderExceeded', ...
                 'Requested evaluation order (%d) exceeds map order (%d). Using map order.', max_eval_order, map_order);
         max_eval_order = map_order;
    elseif max_eval_order < 0
          warning('trace_ray_fan_segmented_da_map:NegativeOrder', ...
                 'Requested evaluation order (%d) is negative. Using order 0.', max_eval_order);
         max_eval_order = 0;
     else
         fprintf('  Using specified evaluation order: %d\n', max_eval_order);
     end
     if S_map_bender{1}.Order ~= map_order || S_steps2{1}{1}.Order ~= map_order
          warning('trace_ray_fan_segmented_da_map:OrderMismatch', ...
                  'DA order mismatch between integration segments or bender map.');
     end
     if abs(z_steps1(end) - z_bend_start) > 1e-6 || abs(z_steps2(1) - z_bend_end) > 1e-6
          warning('trace_ray_fan_segmented_da_map:BoundaryMismatch', ...
                  'Z-boundary mismatch between integration results and bender start/end points.');
     end

    % --- Define Initial Ray Fan Grid ---
    if num_r_rays == 1, r_devs = 0; else r_devs = linspace(-max_r0_dev, max_r0_dev, num_r_rays); end
    if num_s_rays == 1, s_devs = 0; else s_devs = linspace(-max_s0_dev, max_s0_dev, num_s_rays); end
    [R_dev_grid, S_dev_grid] = meshgrid(r_devs, s_devs);
    initial_r0s = center_r0 + R_dev_grid(:);
    initial_s0s = center_s0 + S_dev_grid(:); % s0 = p_r0
    initial_states = [initial_r0s'; initial_s0s']; % 2 x num_rays_total
    initial_states_out = initial_states;
    num_rays_total = size(initial_states, 2);
    fprintf('  Generating initial fan: %d x %d = %d total rays.\n', num_r_rays, num_s_rays, num_rays_total);
    fprintf('    r0 range: [%.3f, %.3f] (Center: %.3f)\n', min(initial_r0s), max(initial_r0s), center_r0);
    fprintf('    p_r0 range: [%.3f, %.3f] (Center: %.3f)\n', min(initial_s0s), max(initial_s0s), center_s0);

    % --- Trace Rays by Evaluating DA Maps Sequentially ---
    num_z_points1 = numel(z_steps1);
    num_z_points2 = numel(z_steps2); % Includes start point z_bend_end
    num_z_points_total = num_z_points1 + num_z_points2; % One point for bender end

    % Preallocate combined storage
    r_positions_combined = zeros(num_z_points_total, num_rays_total);
    z_positions_combined = zeros(num_z_points_total, 1);

    fprintf('  Tracing %d rays by evaluating DA maps (Order <= %d) through %d steps...\n', num_rays_total, max_eval_order, num_z_points_total - 1);
    tic;
    calculation_errors = false;
    X_at_bend_start = zeros(2, num_rays_total); % Store state before bender
    X_at_bend_end = zeros(2, num_rays_total);   % Store state after bender

    % --- Trace Section 1 ---
    fprintf('    Tracing section 1 (%d steps)...\n', num_z_points1);
    for k = 1:num_rays_total
        initial_state_k = initial_states(:, k)'; % Row vector for evaluate
        for i = 1:num_z_points1
            if calculation_errors && any(isnan(r_positions_combined(i:end, k))), continue; end
            S_map = S_steps1{i};
            try
                r_new = S_map{1}.evaluate(initial_state_k, max_eval_order);
                p_new = S_map{2}.evaluate(initial_state_k, max_eval_order); % Also evaluate p
                if ~isfinite(r_new) || ~isfinite(p_new)
                     warning('trace_ray_fan_segmented:EvalNaNInf1', 'NaN/Inf result for ray %d at step %d (z=%.4f). Setting rest to NaN.', k, i, z_steps1(i));
                     r_positions_combined(i:end, k) = NaN;
                     calculation_errors = true; continue;
                 end
                r_positions_combined(i, k) = r_new;
                if i == num_z_points1 % Store state at end of section 1
                    X_at_bend_start(:, k) = [r_new; p_new];
                end
            catch ME_eval
                 fprintf('Error eval DA map (Sec 1) ray %d step %d (z=%.4f): %s\n', k, i, z_steps1(i), ME_eval.message);
                 r_positions_combined(i:end, k) = NaN;
                 calculation_errors = true;
            end
        end % i loop (z steps)
    end % k loop (rays)
    z_positions_combined(1:num_z_points1) = z_steps1;

    % --- Trace Bender Section ---
    fprintf('    Tracing bender section (1 step)...\n');
    bend_step_idx = num_z_points1 + 1; % Index for bender end point
    for k = 1:num_rays_total
         if calculation_errors && any(isnan(r_positions_combined(bend_step_idx:end, k))), continue; end
         state_before_bend_k = X_at_bend_start(:, k)'; % Row vector
         try
             r_new = S_map_bender{1}.evaluate(state_before_bend_k, max_eval_order);
             p_new = S_map_bender{2}.evaluate(state_before_bend_k, max_eval_order);
             if ~isfinite(r_new) || ~isfinite(p_new)
                 warning('trace_ray_fan_segmented:EvalNaNInfBender', 'NaN/Inf result for ray %d at bender step. Setting rest to NaN.', k);
                 r_positions_combined(bend_step_idx:end, k) = NaN;
                 calculation_errors = true; continue;
             end
             r_positions_combined(bend_step_idx, k) = r_new;
             X_at_bend_end(:, k) = [r_new; p_new]; % Store state after bender
         catch ME_eval
             fprintf('Error eval DA map (Bender) ray %d: %s\n', k, ME_eval.message);
             r_positions_combined(bend_step_idx:end, k) = NaN;
             calculation_errors = true;
         end
    end % k loop (rays)
    z_positions_combined(bend_step_idx) = z_bend_end;

    % --- Trace Section 3 ---
    fprintf('    Tracing section 3 (%d steps)...\n', num_z_points2 - 1);
    for k = 1:num_rays_total
        state_after_bend_k = X_at_bend_end(:, k)'; % Row vector for evaluate
        % Loop from the *second* step of integrationResults2, as the first is identity map at z_bend_end
        for i = 2:num_z_points2
            global_idx = bend_step_idx + i - 1; % Index in combined array
            if calculation_errors && any(isnan(r_positions_combined(global_idx:end, k))), continue; end
            S_map = S_steps2{i}; % Cumulative map from z_bend_end
            try
                r_new = S_map{1}.evaluate(state_after_bend_k, max_eval_order);
                p_new = S_map{2}.evaluate(state_after_bend_k, max_eval_order);
                if ~isfinite(r_new) || ~isfinite(p_new)
                     warning('trace_ray_fan_segmented:EvalNaNInf2', 'NaN/Inf result for ray %d at step %d (z=%.4f). Setting rest to NaN.', k, i, z_steps2(i));
                     r_positions_combined(global_idx:end, k) = NaN;
                     calculation_errors = true; continue;
                 end
                r_positions_combined(global_idx, k) = r_new;
            catch ME_eval
                 fprintf('Error eval DA map (Sec 2) ray %d step %d (z=%.4f): %s\n', k, i, z_steps2(i), ME_eval.message);
                 r_positions_combined(global_idx:end, k) = NaN;
                 calculation_errors = true;
            end
        end % i loop (z steps)
    end % k loop (rays)
    z_positions_combined((bend_step_idx+1):end) = z_steps2(2:end); % Fill in z values for section 3

    trace_time = toc;
    fprintf('  Finished ray tracing in %.2f seconds.\n', trace_time);
     if calculation_errors, fprintf('  NOTE: Errors occurred during DA map evaluation.\n'); end

    % --- Plot Ray Fan ---
    figure(2); clf;
    hold on;

    % Color mapping based on initial r0
    initial_r_values = initial_states(1, :);
    min_r0 = min(initial_r_values);
    max_r0 = max(initial_r_values);
    if min_r0 == max_r0, normalized_r0 = ones(1, num_rays_total) * 0.5;
    else, normalized_r0 = (initial_r_values - min_r0) / (max_r0 - min_r0); end
    cmap = cool(256);
    color_indices = round(normalized_r0 * (size(cmap, 1) - 1)) + 1;

    % Plot rays individually
    for k = 1:num_rays_total
        ray_r_positions = r_positions_combined(:, k);
        valid_indices = ~isnan(ray_r_positions);
        if any(valid_indices)
            plot(z_positions_combined(valid_indices), ray_r_positions(valid_indices), ...
                 '-', 'Color', cmap(color_indices(k), :), 'LineWidth', 1.0);
        end
    end

    % Plot electrode geometry (using full stack including bender placeholder geometry if desired)
    if isfield(visParams, 'electrodes_full_stack') && ~isempty(visParams.electrodes_full_stack)
        electrodes_vis = visParams.electrodes_full_stack; % Use the stack with the space
    else
        electrodes_vis = visParams.electrodes; % Fallback to BEM stack
    end
    for k = 1:numel(electrodes_vis)
        elec = electrodes_vis(k);
        % Check if rs and zs are valid before plotting
        if isprop(elec,'rs') && isprop(elec,'zs') && ~isempty(elec.rs) && ~isempty(elec.zs)
            plot(elec.zs, elec.rs, '-k', 'LineWidth', 1.5);
            plot(elec.zs, -elec.rs, '-k', 'LineWidth', 1.5);
        end
    end

    xlabel('Axial Position z');
    ylabel('Radial Position r');
    title(sprintf('Ray Fan (Segmented DA Eval Order <= %d, %dx%d Rays)', max_eval_order, num_r_rays, num_s_rays));
    grid on;
    axis tight;

    ylim_curr = ylim;
    ylim_max_abs = max(abs(ylim_curr));
    ylim_required = max(ylim_max_abs, simParams.R * 1.1);
    ylim([-ylim_required, ylim_required]);

    % Mark boundaries
    xline(z_positions_combined(1), ':k', 'Label', 'Start');
    xline(z_bend_start, '--r', 'Label', 'Bender Start');
    xline(z_bend_end, '--r', 'Label', 'Bender End');
    xline(z_positions_combined(end), ':k', 'Label', 'End');

    % Add colorbar
    colormap(cool);
    cbar = colorbar;
    ylabel(cbar, 'Initial Radial Position (r_0)');
    if min_r0 == max_r0, caxis([min_r0 - 0.1, max_r0 + 0.1]); else, caxis([min_r0, max_r0]); end

    hold off;
    fprintf('  Segmented ray fan plot generated in Figure 2 (colored by initial r0).\n');

end % function trace_ray_fan_segmented_da_map