function [r_positions, initial_states_out] = trace_ray_fan_da_map(traceParams, integrationResults, visParams, max_eval_order)
%trace_ray_fan_da_map Traces a fan of rays using full DA maps evaluated at each step,
%   with flexible fan definition and optional evaluation order truncation.
%   Plots rays color-coded by initial radial position.
%
%   INPUTS:
%       traceParams: Struct with ray tracing parameters:
%                    .num_r_rays (number of rays in position dimension)
%                    .num_s_rays (number of rays in slope dimension)
%                    .max_r0_dev (max deviation from center_r0)
%                    .max_s0_dev (max deviation from center_s0 = p_r = dr/dz)
%                    .center_r0 (central initial radial position)
%                    .center_s0 (central initial radial slope p_r = dr/dz)
%       integrationResults: Struct containing outputs from the integrator:
%                           .z_steps  (vector of z coordinates)
%                           .S_steps  (cell array of DA state maps {r_map, p_map})
%       visParams: Struct with visualization info (electrodes, simParams)
%       max_eval_order: (Optional) Maximum DA order to use for map evaluation.
%                       Defaults to the full order of the stored maps.
%
%   OUTPUTS:
%       r_positions: Matrix (num_z_points x num_rays_total) of radial positions.
%       initial_states_out: Matrix (2 x num_rays_total) of the initial [r0; p0] used.
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

    z_steps = integrationResults.z_steps;
    S_steps = integrationResults.S_steps; % Cell array of {r_map, pr_map} at each z_step

    electrodes = visParams.electrodes;
    simParams = visParams.simParams;

    % --- Input Validation & Defaults ---
    if numel(z_steps) ~= numel(S_steps)
        error('trace_ray_fan_da_map:InputMismatch', 'Mismatch between z_steps and S_steps count.');
    end
    if numel(z_steps) < 2
        error('trace_ray_fan_da_map:NotEnoughSteps', 'Need at least 2 steps for tracing.');
    end
    if ~isa(S_steps{1}{1}, 'DiffAlg') || ~ismethod(S_steps{1}{1}, 'evaluate')
         error('trace_ray_fan_da_map:MissingEvaluate', 'DA objects lack "evaluate" method.');
    end
     map_order = S_steps{1}{1}.Order; % Get DA order from the map itself
     if nargin < 4 || isempty(max_eval_order)
         max_eval_order = map_order; % Default to full map order
         fprintf('  Using default evaluation order: %d (full map order)\n', max_eval_order);
     elseif max_eval_order > map_order
         warning('trace_ray_fan_da_map:OrderExceeded', ...
                 'Requested evaluation order (%d) exceeds map order (%d). Using map order.', max_eval_order, map_order);
         max_eval_order = map_order;
     elseif max_eval_order < 0
          warning('trace_ray_fan_da_map:NegativeOrder', ...
                 'Requested evaluation order (%d) is negative. Using order 0.', max_eval_order);
         max_eval_order = 0;
     else
         fprintf('  Using specified evaluation order: %d\n', max_eval_order);
     end

    % --- Define Initial Ray Fan Grid ---
    if num_r_rays == 1
        r_devs = 0;
    else
        r_devs = linspace(-max_r0_dev, max_r0_dev, num_r_rays);
    end
    if num_s_rays == 1
        s_devs = 0;
    else
        s_devs = linspace(-max_s0_dev, max_s0_dev, num_s_rays);
    end
    [R_dev_grid, S_dev_grid] = meshgrid(r_devs, s_devs);
    initial_r0s = center_r0 + R_dev_grid(:);
    initial_s0s = center_s0 + S_dev_grid(:); % s0 = p_r0
    initial_states = [initial_r0s'; initial_s0s']; % 2 x num_rays_total
    initial_states_out = initial_states;
    num_rays_total = size(initial_states, 2);
    fprintf('  Generating initial fan: %d x %d = %d total rays.\n', num_r_rays, num_s_rays, num_rays_total);
    fprintf('    r0 range: [%.3f, %.3f] (Center: %.3f)\n', min(initial_r0s), max(initial_r0s), center_r0);
    fprintf('    p_r0 range: [%.3f, %.3f] (Center: %.3f)\n', min(initial_s0s), max(initial_s0s), center_s0);

    % --- Trace Rays by Evaluating DA Maps ---
    num_z_points = numel(z_steps);
    r_positions = zeros(num_z_points, num_rays_total);
    p_positions = zeros(num_z_points, num_rays_total); % Optional storage
    r_positions(1, :) = initial_states(1, :);
    p_positions(1, :) = initial_states(2, :);
    fprintf('  Tracing %d rays by evaluating DA maps (Order <= %d) at %d steps...\n', num_rays_total, max_eval_order, num_z_points - 1);
    tic;
    calculation_errors = false;
    for i = 1:(num_z_points - 1)
        S_map = S_steps{i+1};
        r_map = S_map{1};
        p_map = S_map{2};
        for k = 1:num_rays_total
             if calculation_errors && any(isnan(r_positions(i+1:end, k))), continue; end
            initial_state_k = initial_states(:, k)';
            try
                r_new = r_map.evaluate(initial_state_k, max_eval_order);
                p_new = p_map.evaluate(initial_state_k, max_eval_order);
                 if ~isfinite(r_new) || ~isfinite(p_new)
                     warning('trace_ray_fan_da_map:EvalNaNInf', 'NaN/Inf result for ray %d at step %d. Setting rest to NaN.', k, i);
                     r_positions((i+1):end, k) = NaN; p_positions((i+1):end, k) = NaN;
                     calculation_errors = true; continue;
                 end
                r_positions(i+1, k) = r_new;
                p_positions(i+1, k) = p_new;
            catch ME_eval
                 fprintf('Error eval DA map ray %d step %d (z=%.4f): %s\n', k, i, z_steps(i+1), ME_eval.message);
                 r_positions((i+1):end, k) = NaN; p_positions((i+1):end, k) = NaN;
                 calculation_errors = true;
            end
        end % k loop
    end % i loop
    trace_time = toc;
    fprintf('  Finished ray tracing in %.2f seconds.\n', trace_time);
     if calculation_errors, fprintf('  NOTE: Errors occurred during DA map evaluation.\n'); end

    % --- Plot Ray Fan ---
    figure(2); clf;
    hold on;

    % --- NEW: Color mapping based on initial r0 ---
    initial_r_values = initial_states(1, :); % Get the row vector of initial r0s
    min_r0 = min(initial_r_values);
    max_r0 = max(initial_r_values);
    if min_r0 == max_r0 % Avoid division by zero if all r0 are the same
        normalized_r0 = ones(1, num_rays_total) * 0.5; % Assign middle color
    else
        % Normalize r0 to [0, 1] range
        normalized_r0 = (initial_r_values - min_r0) / (max_r0 - min_r0);
    end
    % Create colormap
    cmap = cool(256); % Use 256 colors for smoother gradient
    % Map normalized r0 to colormap indices
    color_indices = round(normalized_r0 * (size(cmap, 1) - 1)) + 1;

    % --- Plot rays individually with specific colors ---
    for k = 1:num_rays_total
        ray_r_positions = r_positions(:, k);
        valid_indices = ~isnan(ray_r_positions); % Find non-NaN points for this ray
        if any(valid_indices) % Only plot if there are valid points
            plot(z_steps(valid_indices), ray_r_positions(valid_indices), ...
                 '-', 'Color', cmap(color_indices(k), :), 'LineWidth', 1.0);
        end
    end
    % --- END NEW color plotting ---

    % Plot electrode geometry
    for k = 1:numel(electrodes)
        elec = electrodes(k);
        plot(elec.zs, elec.rs, '-k', 'LineWidth', 1.5);
        plot(elec.zs, -elec.rs, '-k', 'LineWidth', 1.5);
    end

    xlabel('Axial Position z');
    ylabel('Radial Position r');
    title(sprintf('Ray Fan (Eval Order <= %d, %dx%d Rays)', max_eval_order, num_r_rays, num_s_rays));
    grid on;
    axis tight;

    ylim_curr = ylim;
    ylim_max_abs = max(abs(ylim_curr));
    ylim_required = max(ylim_max_abs, simParams.R * 1.1);
    ylim([-ylim_required, ylim_required]);

    xline(z_steps(1), ':k', 'Label', 'Start');
    xline(z_steps(end), ':k', 'Label', 'End');

    % Add a colorbar to show the mapping (optional but helpful)
    colormap(cool); % Apply the cool colormap to the figure
    cbar = colorbar;
    ylabel(cbar, 'Initial Radial Position (r_0)');
    % Set colorbar limits based on actual r0 range
    if min_r0 == max_r0
        caxis([min_r0 - 0.1, max_r0 + 0.1]); % Adjust if single r0
    else
        caxis([min_r0, max_r0]);
    end

    hold off;

    fprintf('  Ray fan plot generated in Figure 2 (colored by initial r0).\n');

end