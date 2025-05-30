function [r_positions] = trace_ray_fan(traceParams, rk4Results, visParams)
%trace_ray_fan Traces a fan of rays using linear matrices and plots them.
%
%   INPUTS:
%       traceParams: Struct with ray tracing parameters (num_rays, max_r0, initial_p0)
%       rk4Results: Struct with RK4 results (all_R_matrices, z_boundaries)
%       visParams: Struct with visualization info (electrodes, simParams)
%
%   OUTPUTS:
%       r_positions: Matrix (num_steps+1 x num_rays) of radial positions at boundaries.
%
%   Displays Figure 2 with the ray fan plot.

    % --- Extract Parameters ---
    num_rays    = traceParams.num_rays;
    max_r0      = traceParams.max_r0;
    initial_p0  = traceParams.initial_p0;

    all_R_matrices = rk4Results.all_R_matrices;
    z_boundaries   = rk4Results.z_boundaries;

    electrodes = visParams.electrodes;
    simParams = visParams.simParams; % Contains R for plot limits

    % --- Define Initial Ray Fan ---
    if mod(num_rays, 2) == 0
        warning('trace_ray_fan:num_rays', 'num_rays should ideally be odd to include central ray. Using %d rays.', num_rays);
    end
    initial_r0s = linspace(-max_r0, max_r0, num_rays);
    initial_states = [initial_r0s; repmat(initial_p0, 1, num_rays)]; % 2 x num_rays [r0; p0]

    % --- Trace Rays Through Slices ---
    num_steps = numel(all_R_matrices); % Should be num_total_slices
    if num_steps ~= (numel(z_boundaries) - 1)
        error('Mismatch between number of R matrices (%d) and z boundaries (%d).', num_steps, numel(z_boundaries));
    end

    r_positions = zeros(num_steps + 1, num_rays); % Store r at each boundary
    p_positions = zeros(num_steps + 1, num_rays); % Store p at each boundary (optional, not returned)

    % Set initial conditions at z_boundaries(1)
    r_positions(1, :) = initial_states(1, :);
    p_positions(1, :) = initial_states(2, :);
    current_states = initial_states; % State [r; p] at the start of the current slice

    % Apply R matrix for each slice to get state at the end of the slice
    for i = 1:num_steps
        R_slice = all_R_matrices{i};
        current_states = R_slice * current_states; % Apply map: X_end = R_slice * X_start
        r_positions(i+1, :) = current_states(1, :); % Store r at z_boundaries(i+1)
        p_positions(i+1, :) = current_states(2, :); % Store p at z_boundaries(i+1)
    end

    % --- Plot Ray Fan ---
    figure(2); clf; % Create/clear figure 2
    hold on;

    % Plot the traced rays
    plot(z_boundaries, r_positions, '-b'); % Plot r vs z for all rays

    % Plot electrode geometry for context
    for k = 1:numel(electrodes)
        elec = electrodes(k);
        plot(elec.zs, elec.rs, '-k', 'LineWidth', 1.5);  % Upper boundary
        plot(elec.zs, -elec.rs, '-k', 'LineWidth', 1.5); % Lower boundary
    end

    % Add labels and title
    xlabel('Axial Position z');
    ylabel('Radial Position r');
    title(sprintf('Ray Fan Plot (Linear Tracing, %d Rays, p_0=%.2f)', num_rays, initial_p0));
    grid on;
    axis tight; % Adjust axis limits initially to data

    % Set symmetric y-limits, ensuring they encompass the electrodes
    ylim_curr = ylim; % Get current y-limits
    ylim_max_abs = max(abs(ylim_curr)); % Find max absolute y from data
    ylim_required = max(ylim_max_abs, simParams.R * 1.1); % Ensure at least +/- 1.1*R
    ylim([-ylim_required, ylim_required]);

    % Add vertical lines for integration start/end
    xline(z_boundaries(1), '--k', 'Label', 'Start');
    xline(z_boundaries(end), '--k', 'Label', 'End');

    hold off;

    fprintf('  Ray fan plot generated in Figure 2.\n');

end
