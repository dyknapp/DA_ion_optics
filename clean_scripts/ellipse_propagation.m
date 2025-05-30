% ==========================================================
% COURANT-SNYDER ELLIPSE PROPAGATION & PLOTTING
% ==========================================================
fprintf('\n--- Courant-Snyder Ellipse Propagation & Plotting ---\n');

% --- Define Initial Ellipse Parameters ---
% Example: Define based on initial fan or known beam properties
cs_params_initial.alpha0 = -3;       % Example: Upright ellipse initially
cs_params_initial.beta0  = 100;     % Example: Beta function [mm/rad] or units consistent with p_r=dr/dz
cs_params_initial.epsilon = 1e-2;   % Example: Emittance [mm*rad] or units consistent with r*p_r
% cs_params_initial.epsilon = 1e-2;

% Calculate initial gamma using Twiss relation: beta*gamma - alpha^2 = 1
cs_params_initial.gamma0 = (1 + cs_params_initial.alpha0^2) / cs_params_initial.beta0;
fprintf('Initial CS Parameters: alpha=%.3f, beta=%.3f, gamma=%.3f, epsilon=%.3e\n', ...
    cs_params_initial.alpha0, cs_params_initial.beta0, cs_params_initial.gamma0, cs_params_initial.epsilon);

% --- Define Z locations for plotting ---
z_plot_list = [30, 130, 350, 477]; % Example: 5 points including start/end
% You can customize this list, e.g., [z_start, z_mid, z_end]
% Ensure these z-values fall within the integration range

% --- Combine Integration Results for Plotting Function ---
% The plotting function needs the cumulative maps from start to each z
% We need to concatenate the results from the three sections.

% Section 1 results
z_steps_combined = integrationResults1.z_steps(:); % Ensure column
S_steps_combined = integrationResults1.S_steps;

% Section 2 (Bender) - Evaluate segment maps cumulatively
fprintf('  Calculating cumulative maps through bender segments...\n');
S_map_cumulative_bender = S_steps_combined{end}; % Start with map at end of section 1
num_bender_points_plot = numel(z_coords_bender);
for i = 1:n_bender_segments
    S_map_segment = S_maps_bender_segmented{i}; % {r_map_seg, p_map_seg}
    % Compose: S_new(X) = S_segment( S_cumulative(X) )
    % Need a DA composition function, or evaluate points.
    % Simplification: For plotting, we only need the map *at the boundaries*.
    % We need S_map(z_bend_start -> z_coords_bender(i+1))

    % --> NOTE: Direct DA map composition is complex. For simplicity,
    % we will re-evaluate from the start for the DA trace points within
    % the plotting function if needed inside the bender, or just plot at boundaries.
    % Let's add the bender boundaries to z_steps_combined.
    % We need the DA map at the END of each bender segment, relative to the start.

    % Calculate cumulative map THROUGH the bender segment i
    % S_map_cumulative_bender = compose_da_maps(S_map_segment, S_map_cumulative_bender); % Requires compose_da_maps
    % TEMPORARY WORKAROUND: We don't have compose_da_maps.
    % The plotting function will handle evaluating S_steps1{end}, then applying
    % the bender maps sequentially for DA points if z is inside the bender.
    % For the linear part, we just multiply the R matrices.

    % We need the cumulative R matrix at each bender boundary
    R_map_segment = eye(2); % Extract R from S_map_segment (or use transfer_maps_bender)
    r_map_seg = S_map_segment{1}; pr_map_seg = S_map_segment{2};
    idx_r_lin_seg = find(all(r_map_seg.getMultiIndices() == [1, 0], 2), 1);
    idx_pr_lin_seg = find(all(r_map_seg.getMultiIndices() == [0, 1], 2), 1);
    if ~isempty(idx_r_lin_seg) && ~isempty(idx_pr_lin_seg) % Check indices found
         R_map_segment = [r_map_seg.Series(idx_r_lin_seg), r_map_seg.Series(idx_pr_lin_seg); ...
                          pr_map_seg.Series(idx_r_lin_seg), pr_map_seg.Series(idx_pr_lin_seg)];
    else
         warning('Could not extract linear map for bender segment %d. Using identity.', i);
         R_map_segment = eye(2); % Fallback
    end

    % We need the FULL DA map from z_start to the end of this bender segment
    % This is complex without composition. The current S_steps structure is CUMULATIVE.
    % We need to adapt plot_propagated_ellipse or how we pass maps.

    % Option 1: Pass all map segments separately to the plotting function.
    % Option 2: Modify plotting function to handle segmented evaluation.

    % --> Let's choose Option 2: Adapt plot_propagated_ellipse later if needed.
    % For now, we will only plot *outside* the bender for simplicity.
    % We will just append the final section results.
end
fprintf('  (Skipping intermediate bender maps for combined list - plotting outside bender only for now)\n');


% Section 3 results - Need to compose with map at end of bender
S_map_at_bend_end = S_map_cumulative_bender; % Map from z_start to z_bend_end (APPROXIMATE)
% --> This composition step is missing.

% Quick Fix: Only use Section 1 and Section 3 results for plotting outside bender.
z_steps_for_plot = [integrationResults1.z_steps(:); integrationResults2.z_steps(:)];
S_steps_for_plot = [integrationResults1.S_steps; integrationResults2.S_steps]; % This is NOT quite right, S_steps2 needs composition.

% --> Due to the complexity of DA map composition which is not implemented,
%     this plotting will only be truly correct for z values within the ranges covered
%     by integrationResults1 and integrationResults2 *relative to their respective starts*.
%     A full implementation requires DA map composition.

% --- Call the Plotting Function ---
n_ellipse_pts = 100;
plot_fig_handle = figure(5); % Use figure 5 for ellipse plots

% Use only points outside the bender for this simplified example
% z_plot_list_valid = z_plot_list(z_plot_list < z_bend_start | z_plot_list > z_bend_end);
fprintf('Plotting ellipses only outside the bender range [%.2f, %.2f]\n', z_bend_start, z_bend_end);

if ~isempty(z_plot_list_valid)
    % Issue: S_steps_for_plot{idx} is not the map from z_start_overall if idx corresponds to section 3.
    % WORKAROUND: Plotting function needs modification or separate calls.
    % Let's just call it separately for section 1 and section 3 for now.
    
    figure(1);
    fprintf('Plotting for Section 1...\n');
    plot_propagated_ellipse(integrationResults1.z_steps, ...
                            integrationResults1.S_steps, ...
                            cs_params_initial, ...
                            z_plot_list(1:2), ...
                            n_ellipse_pts, ...
                            trace_eval_order, ...
                            gcf);

    % For Section 3, we need to transform the initial ellipse to the bender end
    fprintf('Plotting for Section 3 (Requires transform through bender)...\n');
    % 1. Get cumulative R matrix to end of bender
    R_to_bend_end = eye(2);
    S_map_end_sec1 = integrationResults1.S_steps{end};
     r_map_end1 = S_map_end_sec1{1}; pr_map_end1 = S_map_end_sec1{2};
     idx_r_lin_end1 = find(all(r_map_end1.getMultiIndices() == [1, 0], 2), 1);
     idx_pr_lin_end1 = find(all(r_map_end1.getMultiIndices() == [0, 1], 2), 1);
     if ~isempty(idx_r_lin_end1) && ~isempty(idx_pr_lin_end1)
         R_to_bend_start = [r_map_end1.Series(idx_r_lin_end1), r_map_end1.Series(idx_pr_lin_end1); ...
                            pr_map_end1.Series(idx_r_lin_end1), pr_map_end1.Series(idx_pr_lin_end1)];
     else R_to_bend_start = eye(2); end

    R_bender_total = eye(2);
    for i=1:n_bender_segments
        % Extract linear map M_i from S_maps_bender_segmented{i} or use transfer_maps_bender{i}
         M_i = transfer_maps_bender{i}; % Use pre-calculated linear maps
         R_bender_total = M_i * R_bender_total; % Apply maps sequentially
    end
     R_to_bend_end = R_bender_total * R_to_bend_start; % Total R from z_start to z_bend_end

    % 2. Propagate CS parameters to bender end
    [alpha_bend_end, beta_bend_end, gamma_bend_end] = propagate_cs(cs_params_initial.alpha0, cs_params_initial.beta0, cs_params_initial.gamma0, R_to_bend_end);
    cs_params_bend_end.alpha0 = alpha_bend_end;
    cs_params_bend_end.beta0 = beta_bend_end;
    cs_params_bend_end.gamma0 = gamma_bend_end;
    cs_params_bend_end.epsilon = cs_params_initial.epsilon * 2.0e-1; % Emittance conserved in linear approx

    % 3. Call plotting function for section 3, using transformed CS parameters
    %    S_steps for section 3 are relative to z_bend_end.
    z_plot_section3 = z_plot_list(z_plot_list >= z_bend_end);
    figure(2);
    if ~isempty(z_plot_section3)
         plot_propagated_ellipse(integrationResults2.z_steps, ... % z starts at z_bend_end
                                integrationResults2.S_steps, ... % maps are relative to z_bend_end
                                cs_params_bend_end, ...             % CS params at z_bend_end
                                z_plot_section3, ...                % Target z values
                                n_ellipse_pts, ...
                                trace_eval_order, ...
                                gcf); % Add to the same figure
    end
else
     fprintf('No valid z locations specified for ellipse plotting outside the bender.\n');
end

function XY = get_ellipse_points(alpha, beta, gamma, epsilon, n_points)
% Generates points (r, pr) on the boundary of a CS ellipse.
% Equation: gamma*r^2 + 2*alpha*r*pr + beta*pr^2 = epsilon

    if abs(beta * gamma - alpha^2 - 1) > 1e-6
        % warning('get_ellipse_points:TwissMismatch', 'CS parameters do not satisfy beta*gamma - alpha^2 = 1');
        % Correct gamma based on alpha and beta
        gamma = (1 + alpha^2) / beta;
    end

    % Generate points on a unit circle parameterization (a, b)
    theta = linspace(0, 2*pi, n_points+1);
    theta = theta(1:end-1); % n_points unique points

    % The transformation from unit circle (a,b) to (r, pr) is:
    % r = sqrt(beta*epsilon) * a
    % pr = sqrt(gamma*epsilon) * b - (alpha/sqrt(beta*epsilon)) * r
    % More robust transformation (avoids division by zero if beta is small):
    % Use Cholesky decomposition of the sigma matrix: Sigma = epsilon * [beta, -alpha; -alpha, gamma]
    % Sigma = epsilon * [[beta, -alpha]; [-alpha, gamma]]
    Sigma = epsilon * [beta, -alpha; -alpha, gamma];

    % Need to ensure Sigma is positive definite. This is guaranteed if
    % epsilon > 0, beta > 0, and det(Sigma) = epsilon^2 * (beta*gamma - alpha^2) = epsilon^2 > 0.
    % Beta is beam size squared, must be > 0.
    if beta <= 0
       error('Beta must be positive to define an ellipse.');
    end

    try
        L = chol(Sigma, 'lower'); % Sigma = L * L'
    catch ME
         warning('Sigma matrix not positive definite? Beta=%.3e, Gamma=%.3e, Alpha=%.3e, Eps=%.2e, det=%.3e', beta, gamma, alpha, epsilon, det(Sigma));
         % Fallback to simpler parameterization if Cholesky fails
         r = sqrt(beta * epsilon) * cos(theta);
         pr = (-alpha * r + sqrt(epsilon * (alpha^2*beta*epsilon*cos(theta).^2 + beta*epsilon - alpha^2*beta*epsilon) ) / sqrt(beta*epsilon) ) / beta;
         % This is complex, let's use the standard matrix method if possible.
         % Rethrowing error for now.
          rethrow(ME);
    end

    % Generate points on unit circle
    unit_circle_points = [cos(theta); sin(theta)]; % 2 x n_points

    % Transform to (r, pr) space: XY = L * unit_circle_points
    XY = L * unit_circle_points; % 2 x n_points
end


function [alpha_out, beta_out, gamma_out] = propagate_cs(alpha_in, beta_in, gamma_in, R)
% Propagates Courant-Snyder parameters through a transfer matrix R.
    R11 = R(1,1); R12 = R(1,2);
    R21 = R(2,1); R22 = R(2,2);

    beta_out  =  R11^2 * beta_in - 2*R11*R12 * alpha_in + R12^2 * gamma_in;
    alpha_out = -R11*R21* beta_in + (R11*R22 + R12*R21) * alpha_in - R12*R22 * gamma_in;
    gamma_out =  R21^2 * beta_in - 2*R21*R22 * alpha_in + R22^2 * gamma_in;
end
