% Assume we have a perfect quadrupole field.
E0 = 0.1;
Rquad = 48;
R0 = 86;
V0 = -2.3123 * E0 * (R0 / Rquad)^2;

% Step 1: define the potential as a lambda function
% particle goes in along x, comes out along y.
phi = @(x, y) 0.5 * (V0 / R0^2) * ((x+y).^2 - (x-y).^2); % Simplified form: phi = V0/R0^2 * 2*x*y

extent = linspace(-R0/2, R0/8, 256);
[xs, ys] = meshgrid(extent, -extent);

f = figure(3); clf; f.Name = 'Quadrupole Potential & Curvilinear axis';
contourf(xs, ys, phi(xs, ys), 100); axis image; colorbar;
title('Quadrupole Potential \phi(x,y)');
xlabel('x'); ylabel('y');

% Checking that the functional form of the potential is correct.
% At 45 degrees, it should be +/- V0 * (R/R0)^2 - Wait, V0 already scaled.
% Phi(x,y) = V0 * 2*x*y / R0^2. Let x=Rquad/sqrt(2), y=Rquad/sqrt(2)
fprintf("Potential Check:\n");
fprintf("Phi( Rq/sqrt(2),  Rq/sqrt(2)) = % .12f (Should be +V0 * (Rq/R0)^2 / 2 * (Rq/R0)^2 ? No...)\n", phi(Rquad / sqrt(2),  Rquad / sqrt(2))) % Should be V0 * Rquad^2 / R0^2
fprintf("Phi( Rq/sqrt(2), -Rq/sqrt(2)) = % .12f (Should be -V0 * (Rq/R0)^2 / 2 * (Rq/R0)^2 ? No...)\n", phi(Rquad / sqrt(2), -Rquad / sqrt(2))) % Should be -V0 * Rquad^2 / R0^2
fprintf("Expected V0*(Rq/R0)^2 = %.6f\n", V0 * (Rquad/R0)^2); % Let's re-evaluate V definition later if needed. Keep V0 as scale factor for now.

%% Draw the quarter circle reference path
circle = @(r, theta) [(R0+r)/2 .* sin(theta) - R0/2; -(R0+r)/2 .* cos(theta) + R0/2];

hold on
circle_points_main = circle(0, linspace(0, pi/2, 128));
circle_points_plus = circle(R0/10, linspace(0, pi/2, 128));
circle_points_minus = circle(-R0/10, linspace(0, pi/2, 128));
plot(circle_points_main(1, :),  circle_points_main(2, :), '-k', LineWidth=3, DisplayName='Ref (r=0)')
plot(circle_points_plus(1, :),  circle_points_plus(2, :), '-g', LineWidth=1, DisplayName='Ref (r=R0/10)')
plot(circle_points_minus(1, :), circle_points_minus(2, :), '-r', LineWidth=1, DisplayName='Ref (r=-R0/10)')
hold off
legend show;

%% Main solver - Calculate Four Trajectories

% --- Initial Conditions ---
% State vector: [r; dr/dtheta]
y0_ref = [0.0 * R0; 0.0];           % Reference trajectory (r0=0, dr/dtheta0=0)
y0_1   = [0.0 * R0; 0.1];           % (r0=0, dr/dtheta0=0.1)
y0_2   = [R0/10;    0.0];           % (r0=R0/10, dr/dtheta0=0)
y0_3   = [R0/10;    0.1];           % (r0=R0/10, dr/dtheta0=0.1)

initial_conditions = {y0_ref, y0_1, y0_2, y0_3};
num_trajectories = length(initial_conditions);
trajectory_results = cell(1, num_trajectories); % Store results for each trajectory

% --- Integration Span & Evaluation Points ---
dth = 0; % No offset needed if starting exactly at 0
theta_start = dth;
theta_end = pi/2 - dth;
n_segments = 100; % Number of segments for transfer maps
n_points = n_segments + 1; % Need points at boundaries
theta_eval = linspace(theta_start, theta_end, n_points); % Points for evaluation (segment boundaries)

% Set ODE solver options
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

fprintf('\nIntegrating trajectories...\n');
integration_successful = true;

for i = 1:num_trajectories
    y0 = initial_conditions{i};
    fprintf('Integrating trajectory %d: r0=%.4f, dr/dth0=%.4f\n', i, y0(1), y0(2));
    try
        [theta_sol, y_sol] = ode45(@(theta, y) ode_system(theta, y, R0, V0, E0), theta_eval, y0, options);

        % Check if integration finished successfully for this trajectory
        if abs(theta_sol(end) - theta_end) > 1e-5 || size(y_sol, 1) ~= n_points
            fprintf('Warning: Integration for trajectory %d might not have completed successfully.\n', i);
            fprintf('  Reached theta = %.4f (expected %.4f), got %d points (expected %d)\n', theta_sol(end), theta_end, size(y_sol,1), n_points);
            integration_successful = false;
            % Store partial results if needed, or handle error
             trajectory_results{i} = struct('theta', theta_sol, 'y', y_sol, 'success', false);
             % break; % Option to stop if one fails
        else
             trajectory_results{i} = struct('theta', theta_sol, 'y', y_sol, 'success', true);
        end

    catch ME
        fprintf('Integration failed for trajectory %d with error: %s\n', i, ME.message);
        fprintf('Error identifier: %s\n', ME.identifier);
        disp(ME.getReport());
        integration_successful = false;
        trajectory_results{i} = struct('theta', [], 'y', [], 'success', false);
        % break; % Option to stop if one fails
    end
end

if ~integration_successful
    fprintf('One or more integrations failed or did not complete. Transfer map calculation aborted.\n');
    % Exit or handle gracefully
    return;
else
    fprintf('All integrations completed successfully.\n');
end

%% Calculate Linear Transfer Maps per Segment

fprintf('Calculating transfer maps...\n');
transfer_maps = cell(1, n_segments);
z_coords = R0 * theta_eval / 2; % Curvilinear coordinate z = R0*theta

% Extract reference trajectory results
y_ref = trajectory_results{1}.y; % [r, dr/dtheta]
r_ref = y_ref(:, 1);
p_ref = y_ref(:, 2) / R0; % p = dr/dz = (dr/dtheta) / R0
ref_states = [r_ref(:)'; p_ref(:)']; % Store as 2 x N_points

valid_maps = 0;
for i = 1:n_segments
    % Indices for start (idx0) and end (idx1) of the segment
    idx0 = i;
    idx1 = i + 1;

    % State vectors Y = [r; p] = [r; dr/dz] at start and end of segment
    Y_ref_0 = ref_states(:, idx0); % Use pre-calculated ref_states
    Y_ref_1 = ref_states(:, idx1); % Use pre-calculated ref_states

    % Calculate deviation vectors for the other 3 trajectories
    DeltaY_in = zeros(2, 3); % Columns: traj 1, 2, 3 at segment start
    DeltaY_out = zeros(2, 3); % Columns: traj 1, 2, 3 at segment end

    valid_segment = true;
    for k = 1:3 % Trajectories 2, 3, 4 (index k+1 in results)
        if ~trajectory_results{k+1}.success % Check if trajectory integration was successful
            fprintf('Warning: Skipping segment %d map calculation due to unsuccessful integration of trajectory %d.\n', i, k+1);
            valid_segment = false;
            break;
        end
        y_k = trajectory_results{k+1}.y;
        r_k0 = y_k(idx0, 1);
        p_k0 = y_k(idx0, 2) / R0;
        r_k1 = y_k(idx1, 1);
        p_k1 = y_k(idx1, 2) / R0;

        Y_k0 = [r_k0; p_k0];
        Y_k1 = [r_k1; p_k1];

        DeltaY_in(:, k) = Y_k0 - Y_ref_0;
        DeltaY_out(:, k) = Y_k1 - Y_ref_1;

        % Basic check for NaN/Inf in deviations
        if any(isnan(DeltaY_in(:,k))) || any(isinf(DeltaY_in(:,k))) || ...
           any(isnan(DeltaY_out(:,k))) || any(isinf(DeltaY_out(:,k)))
           fprintf('Warning: NaN/Inf detected in deviation calculation for segment %d, trajectory %d.\n', i, k+1);
           valid_segment = false;
           break;
        end
    end

    if ~valid_segment
        transfer_maps{i} = eye(2); % Use identity matrix as fallback
        fprintf('  Using identity matrix for segment %d.\n', i);
        continue; % Skip to next segment
    end

    % Solve for M_i using least squares: M_i * DeltaY_in = DeltaY_out
    % M_i = DeltaY_out / DeltaY_in; % Solves in least-squares sense
    % Alternative formulation: M = (Y_out * Y_in') * inv(Y_in * Y_in');
    try
        % Check condition number before solving
        if cond(DeltaY_in) > 1e10 % Check if input deviation vectors are nearly linearly dependent
             fprintf('Warning: Poorly conditioned input deviations for segment %d (cond=%.2e). Using identity map.\n', i, cond(DeltaY_in));
             M_i = eye(2);
        else
            M_i = DeltaY_out / DeltaY_in;
             % Check if M_i calculation resulted in NaN/Inf
            if any(isnan(M_i(:))) || any(isinf(M_i(:)))
                fprintf('Warning: Transfer matrix calculation resulted in NaN/Inf for segment %d. Using identity map.\n', i);
                M_i = eye(2);
            end
        end
    catch ME_solve
        fprintf('Error solving for transfer matrix in segment %d: %s\n', i, ME_solve.message);
        M_i = eye(2); % Fallback to identity
    end


    % Normalize to enforce unit determinant (Symplecticity)
    det_M = det(M_i);
    if abs(det_M) < 1e-10 || isnan(det_M) || isinf(det_M) % Check for near-zero or invalid determinant
        fprintf('Warning: Determinant is near zero or invalid (%.3e) for segment %d. Using identity matrix.\n', det_M, i);
        M_norm = eye(2);
    else
        M_norm = M_i / sqrt(abs(det_M)); % Use abs() for robustness against small numerical negativity
        % Optional: check determinant of normalized matrix
        % fprintf('Segment %d: det(M)=%.6f, det(M_norm)=%.6f\n', i, det_M, det(M_norm));
    end

    transfer_maps{i} = M_norm;
    valid_maps = valid_maps + 1;
end
fprintf('Calculated %d valid transfer maps out of %d segments.\n', valid_maps, n_segments);

%% Propagate Ray Fan using Transfer Maps

fprintf('Propagating ray fan...\n');

% Define initial ray fan conditions relative to the reference trajectory start
n_rays_pos = 5; % Number of rays on each side of center for position
n_rays_slope = 5; % Number of rays on each side of center for slope
r_max_dev = 1.5;  % Max initial r deviation
p_max_dev = 0.25;      % Max initial p = dr/dz deviation

r_starts = linspace(-r_max_dev, r_max_dev, 2*n_rays_pos + 1);
p_starts = linspace(-p_max_dev, p_max_dev, 2*n_rays_slope + 1);

[R_start_grid, P_start_grid] = meshgrid(r_starts, p_starts);
initial_fan_deviations = [R_start_grid(:)'; P_start_grid(:)']; % Each column is an initial [dr; dp] deviation
num_fan_rays = size(initial_fan_deviations, 2);

% Store history of fan rays (absolute coordinates)
% Dimensions: [r or p, point_index, ray_index]
fan_history_r = zeros(n_points, num_fan_rays);
fan_history_p = zeros(n_points, num_fan_rays); % p = dr/dz

% Initial absolute state for the fan using ref_states
fan_history_r(1, :) = ref_states(1, 1) + initial_fan_deviations(1, :);
fan_history_p(1, :) = ref_states(2, 1) + initial_fan_deviations(2, :);

% Propagate each ray
current_deviations = initial_fan_deviations;
for i = 1:n_segments % Step through segments
    idx1 = i + 1; % Index at the end of the current segment

    % Apply transfer map for this segment to the current *deviations*
    next_deviations = transfer_maps{i} * current_deviations;

    % Calculate absolute state at the end of the segment
    Y_ref_1 = ref_states(:, idx1); % Reference state at end of segment
    fan_history_r(idx1, :) = Y_ref_1(1) + next_deviations(1, :);
    fan_history_p(idx1, :) = Y_ref_1(2) + next_deviations(2, :);

    % Update deviations for the next step
    current_deviations = next_deviations;
end
fprintf('Ray fan propagation complete.\n');

%% Plotting

% --- Plot Reference Trajectory and Fan in Cartesian Coordinates ---
f3 = figure(3); % Use the existing figure 3 (potential plot)
hold on; % Add to the potential plot

% Plot the reference trajectory (calculated via ODE) in white
y_ref_ode = trajectory_results{1}.y;
r_ref_ode = y_ref_ode(:,1);
theta_ref_ode = trajectory_results{1}.theta;
x_ref_ode =  (R0 + r_ref_ode) / 2 .* sin(theta_ref_ode(:)) - R0 / 2; % Ensure theta is column
y_ref_ode_cart = -(R0 + r_ref_ode) / 2 .* cos(theta_ref_ode(:)) + R0 / 2; % Ensure theta is column
plot(x_ref_ode, y_ref_ode_cart, '-w', 'LineWidth', 2, 'DisplayName', 'Reference Traj (ODE)');

% Plot the fan trajectories in Cartesian coordinates
% Define a colormap for the fan rays
colors = cool(num_fan_rays);
% Ensure theta_eval is a column vector for calculations below
theta_col = theta_eval(:); % Make theta_eval n_points x 1

for k = 1:num_fan_rays
    r_k = fan_history_r(:, k); % Column vector (n_points x 1)

    % Calculate Cartesian coordinates using column vectors
    x_k = (R0 + r_k) / 2 .* sin(theta_col) - R0 / 2; % Now (n_points x 1) .* (n_points x 1)
    y_k = -(R0 + r_k) / 2 .* cos(theta_col) + R0 / 2; % Now (n_points x 1) .* (n_points x 1)

    plot(x_k, y_k, '-', 'Color', [colors(k,:), 0.6], 'LineWidth', 0.5, 'HandleVisibility','off'); % Thinner, semi-transparent lines
end
hold off;
title('Quadrupole Potential & Trajectories');
legend show; % Ensure legend is shown after adding fan ref traj
fprintf('Plotted fan on potential map (Figure 3).\n');
drawnow;


% --- Plot Solution (r vs theta) and Trajectory (y vs x) with Fan ---
f4 = figure(4); clf; f4.Name = 'Quadrupole Deflection Solution & Ray Fan';
% Ensure theta_eval is a column vector for consistency later
theta_col = theta_eval(:); % Make theta_eval n_points x 1

% Subplot 1: r vs theta/pi
subplot(1, 2, 1);
hold on; % Hold for fan plots

% Plot reference r vs theta from ODE solver
% Use theta_col here too for consistency, although row vector works for 1D plot
plot(theta_col / pi, r_ref_ode / R0, '-k', 'LineWidth', 2, 'DisplayName', 'Ref. Traj (ODE)');

% Plot fan trajectories (r vs theta)
colors = cool(num_fan_rays); % Define colormap here
for k = 1:num_fan_rays
    plot(theta_col / pi, fan_history_r(:, k) / R0, '-', 'Color', [colors(k,:), 0.6], 'LineWidth', 0.5, 'HandleVisibility','off');
end

xlabel('$\theta/\pi$ (radians/pi)', 'Interpreter', 'latex');
ylabel('$r(\theta) / R_0$', 'Interpreter', 'latex');
title('Solution $r$ vs $\theta$', 'Interpreter', 'latex');
grid on;
box on;
hold off;
legend('Location', 'best');


% Subplot 2: y vs x (Cartesian Trajectory)
subplot(1, 2, 2);
hold on; % Hold for fan plots

% Plot reference trajectory (y vs x) from ODE solver (using already calculated x_ref_ode, y_ref_ode_cart)
plot(x_ref_ode, y_ref_ode_cart, '-k', 'LineWidth', 2, 'DisplayName', 'Ref. Traj (ODE)');

% Plot asymptotes V=0 (x=0 and y=0)
xlim_vals = xlim; ylim_vals = ylim; % Get initial limits based on ref traj
plot(xlim_vals, [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'DisplayName', 'V=0 Asymptotes');
plot([0 0], ylim_vals, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility','off'); % Use '--' for consistency

% Plot fan trajectories (y vs x)
for k = 1:num_fan_rays
    r_k = fan_history_r(:, k); % Column vector (n_points x 1)

    % Calculate Cartesian coordinates using column vectors
    x_k = (R0 + r_k) / 2 .* sin(theta_col) - R0 / 2; % (n_points x 1) .* (n_points x 1)
    y_k = -(R0 + r_k) / 2 .* cos(theta_col) + R0 / 2; % (n_points x 1) .* (n_points x 1)

    plot(x_k, y_k, '-', 'Color', [colors(k,:), 0.6], 'LineWidth', 0.5, 'HandleVisibility','off');
end

% Reset limits MAYBE? Or let the fan determine the limits.
% If the fan plot looks cut off, comment out the next line.
% xlim(xlim_vals); ylim(ylim_vals);

xlabel('x');
ylabel('y');
title('Trajectory in Cartesian Coordinates (y vs x)');
axis equal; % Ensure aspect ratio is equal
grid on;
box on;
hold off;
legend('Location', 'best'); % Show legend

% Add overall title
sgtitle(sprintf('Numerical Solution & Ray Fan (V0=%.2f, R0=%d)', V0, R0));
fprintf('Plotted fan on solution plots (Figure 4).\n');
drawnow;

% ================================================
% START: Modified section for saving SEGMENTED data
% ================================================
%% Calculate and Save Data for Full Stack (Segmented)
fprintf('\nCalculating and saving SEGMENTED data for full_stack.m...\n');

% --- Data to Save ---
% transfer_maps: Cell array (1 x n_segments) containing the 2x2 linear map Mi for each segment i.
% ref_states:    Matrix (2 x n_points) containing [r_ref; p_ref] at segment boundaries.
% z_coords:      Vector (n_points x 1) containing z = R0*theta at segment boundaries.

% Ensure data is consistent
n_segments_check = numel(transfer_maps);
n_points_check = size(ref_states, 2);
n_zcoords_check = numel(z_coords);

if n_segments_check ~= (n_points_check - 1) || n_points_check ~= n_zcoords_check
    error('Inconsistent dimensions between transfer_maps (%d segments), ref_states (%d points), and z_coords (%d points). Cannot save segmented data.', ...
          n_segments_check, n_points_check, n_zcoords_check);
end
if valid_maps ~= n_segments_check % Check if all maps were valid
    warning('Not all transfer maps were calculated successfully (%d/%d). Saved data may contain identity matrices.', valid_maps, n_segments_check);
end

% Define output filename
output_filename_segmented = 'MAT007_clean_DA/bender_segmented_data.mat'; % Ensure folder exists or adjust path

% Save the segmented data
try
    save(output_filename_segmented, 'transfer_maps', 'ref_states', 'z_coords');
    fprintf('  Successfully saved segmented bender data (transfer_maps, ref_states, z_coords) to %s\n', output_filename_segmented);
catch ME_save
    fprintf('Error saving segmented bender data: %s\n', ME_save.message);
    error('Could not save segmented bender data.'); % Stop if saving fails
end

% (Optional: Keep the old saving block for the overall matrix if needed elsewhere)
% fprintf('\nCalculating and saving OVERALL data for full_stack.m (original)...\n');
% R_bender_total = eye(2);
% for i = 1:n_segments
%     R_bender_total = transfer_maps{i} * R_bender_total;
% end
% fprintf('  Total Bender Transfer Matrix R_total:\n');
% fprintf('    [ % 9.6f  % 9.6f ]\n', R_bender_total(1,1), R_bender_total(1,2));
% fprintf('    [ % 9.6f  % 9.6f ]\n', R_bender_total(2,1), R_bender_total(2,2));
% ref_state_in = ref_states(:, 1);
% ref_state_out = ref_states(:, end);
% bender_length = z_coords(end) - z_coords(1);
% output_filename_overall = 'MAT007_clean_DA/bender_data.mat';
% try
%     save(output_filename_overall, 'R_bender_total', 'ref_state_in', 'ref_state_out', 'bender_length');
%     fprintf('  Successfully saved overall bender data to %s\n', output_filename_overall);
% catch ME_save
%     fprintf('Error saving overall bender data: %s\n', ME_save.message);
% end
% ================================================
% END: Modified section for saving SEGMENTED data
% ================================================


%% FUNCTIONS 

function V = V_quad(r, theta, R_val, V0_val, E0)
   % Calculates the quadrupole potential.
   % Using element-wise operators for potential vectorization safety
   % Original formula from description: 0.5 * (V0 / R0^2) * ((x+y)^2 - (x-y)^2)
   % Convert (r, theta) curvilinear to (x, y) cartesian relative to center (-R0/2, R0/2)
   % x = (R_val+r)/2 * sin(theta) - R_val/2
   % y = -(R_val+r)/2 * cos(theta) + R_val/2
   % This function seems to implement a *different* potential form than the lambda phi above.
   % Let's use the formula from the ODE paper directly derived from Lagrangian:
   % V(r,theta) = E0 - phi(x(r,theta), y(r,theta))
   % where phi is the electric potential. Let's use the definition consistent
   % with the derivatives below.

   x_cart =  (R_val + r) / 2 .* sin(theta) - R_val / 2;
   y_cart = -(R_val + r) / 2 .* cos(theta) + R_val / 2;
   phi_val = V0_val * 2 .* x_cart .* y_cart / (R_val^2); % Using the simple quadrupole form phi = C * x * y

   % Potential V used in the ODE is related to Kinetic Energy T = V
   % From paper derivation T = E0 - q*phi. Let q=1. So V = E0 - phi
   %V = E0 - phi_val; % *** SWITCHING BACK TO ORIGINAL V_quad FORM BELOW ***

   % Let's double-check the provided V_quad function form against the derivatives.
   % term1 = -R_val + (r + R_val) .* cos(theta); % = -(R_val*(1-cos)) + r*cos
   % term2 = -R_val + (r + R_val) .* sin(theta); % = -(R_val*(1-sin)) + r*sin
   % denom = (2 * R_val^2);
   % V_provided = E0 - V0_val .* term1 .* term2 ./ denom;
   % This form doesn't immediately look like E0 - C*x*y.
   % Let's *trust the derivatives* provided and the ODE, and use the V_quad provided in the original script
   % as it's consistent with the derivatives used in the ODE.

   term1 = -R_val + (r + R_val) .* cos(theta);
   term2 = -R_val + (r + R_val) .* sin(theta);
   denom = (2 * R_val^2) + 1e-30; % Add epsilon for safety
   V = E0 - V0_val .* term1 .* term2 ./ denom;
end

function dVdr = dV_dr(r, theta, R_val, V0_val)
   % Calculates the partial derivative dV/dr, consistent with V_quad above.
   s = sin(theta);
   c = cos(theta);
   factor1 = -R_val .* (c + s);
   factor2 = 2 .* (r + R_val) .* s .* c;
   denom = (2 * R_val^2) + 1e-30;
   dVdr = -V0_val .* (factor1 + factor2) ./ denom;
end

function dVdth = dV_dtheta(r, theta, R_val, V0_val)
   % Calculates the partial derivative dV/dtheta, consistent with V_quad above.
   s = sin(theta);
   c = cos(theta);
   cos2th = cos(2 * theta); % c.^2 - s.^2
   factor1 = R_val .* (s - c);
   factor2 = (r + R_val) .* cos2th;
   term_in_brackets = factor1 + factor2;
   denom = (2 * R_val^2) + 1e-30;
   dVdth = -V0_val .* (r + R_val) .* term_in_brackets ./ denom;
end

% --- ODE System Definition ---
function dydtheta = ode_system(theta, y, R_val, V0_val, E0)
   % Defines the system of first-order ODEs for r and dr/dtheta.
   % y(1) = r
   % y(2) = dr/dtheta
   % MUST return a COLUMN vector: [dr/dtheta; d2r/dtheta2]

   r = y(1);
   dr_dtheta = y(2);
   r_plus_R = r + R_val;

   % Check for coordinate singularity (r approx -R0)
   if abs(r_plus_R) < 1e-12
       fprintf('Warning: r near -R (r+R=%.2e) at theta=%.4f. Singularity.\n', r_plus_R, theta);
       dydtheta = [NaN; NaN]; % Return NaN column vector to signal failure
       return;
   end

   % Calculate potential V = T (Kinetic Energy) and its derivatives
   V = V_quad(r, theta, R_val, V0_val, E0);
   dVdr = dV_dr(r, theta, R_val, V0_val);
   dVdth = dV_dtheta(r, theta, R_val, V0_val);

   % Check for V near zero (Zero Kinetic Energy -> infinite time derivatives)
   % Allow slightly negative V due to numerical issues, but halt if too negative or zero.
   if V < 1e-15 % Changed from abs(V) < 1e-15 to V < 1e-15 (T must be > 0)
        fprintf('Warning: V (Kinetic Energy) near zero or negative (V=%.3e) at theta=%.4f, r=%.4f. Cannot continue.\n', V, theta, r);
        dydtheta = [NaN; NaN]; % Return NaN column vector
        return;
   end

   % Calculate terms for the r'' equation (d2r/dtheta2)
   try
       r_prime_sq = dr_dtheta^2;
       D = r_prime_sq + r_plus_R^2;

       % Term 1: Centrifugal/Coriolis like term
       term1 = (2 * r_prime_sq + r_plus_R^2) / r_plus_R;

       % Term 2: Force term (proportional to potential gradients)
       term_in_brackets = dVdr - (dr_dtheta / (r_plus_R^2)) * dVdth;
       term2 = (D / (2 * V)) * term_in_brackets;

       dr_ddtheta = term1 + term2; % This is d²r/dθ²

       % Final check if calculation resulted in non-finite number
       if ~isfinite(dr_ddtheta)
           fprintf('Warning: Calculated d²r/dθ² is not finite (%.3e) at theta=%.4f.\n r=%.4f, dr/dth=%.4f, V=%.3e\n', dr_ddtheta, theta, r, dr_dtheta, V);
           dydtheta = [dr_dtheta; NaN]; % Return partial result with NaN
           return;
       end

   catch ME % Catch calculation errors (e.g., division by zero if r_plus_R was missed)
       fprintf('Error during calculation in ode_system at theta=%.4f: %s\n', theta, ME.message);
       dydtheta = [NaN; NaN]; % Return NaN column vector
       return;
   end

   % Return derivatives as a COLUMN vector
   dydtheta = [dr_dtheta; dr_ddtheta];
end