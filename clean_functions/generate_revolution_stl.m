function generate_revolution_stl(r_profile, z_profile, thickness, filename, n_revolution)
%GENERATE_REVOLUTION_STL Generates an STL file for a solid of revolution.
%   Generates an STL file representing a solid formed by revolving a 2D
%   profile around the Z-axis. The input profile defines the OUTER surface,
%   and the solid is assumed to have a constant radial thickness.
%
%   Args:
%       r_profile (vector): Vector of radial coordinates (r >= 0) for the
%                           OUTER profile boundary.
%       z_profile (vector): Vector of corresponding axial (z) coordinates
%                           for the OUTER profile boundary. Must be the
%                           same size as r_profile. The pair (r_profile,
%                           z_profile) defines the path to be revolved.
%       thickness (scalar): The constant radial thickness of the solid wall.
%                           Must be positive and generally smaller than min(r_profile).
%       filename (char/string): The name for the output STL file (e.g., 'my_shape.stl').
%       n_revolution (scalar, optional): The number of angular segments used
%                                        for the revolution (controls smoothness).
%                                        Defaults to 50 if empty or not provided.
%
%   Requires:
%       stlwrite function (available from MATLAB File Exchange or built-in
%       in some newer versions/toolboxes).
%
%   Example:
%       % Define a vase-like outer profile
%       z_pts = linspace(0, 10, 20);
%       r_pts = 3 + 2*sin(pi*z_pts/10); % Bulging shape
%       thickness = 0.2;
%       output_file = 'vase.stl';
%       generate_revolution_stl(r_pts, z_pts, thickness, output_file, 60);
%       % You can then view 'vase.stl' in an STL viewer or import into CAD.

    % --- Input Handling & Defaults ---
    if nargin < 5 || isempty(n_revolution)
        n_revolution = 50; % Default number of segments around Z
    end
    if ~isscalar(n_revolution) || n_revolution <= 2
        warning('n_revolution should be a scalar > 2. Setting to default 50.');
        n_revolution = 50;
    end
    n_revolution = round(n_revolution); % Ensure integer

    if ~isvector(r_profile) || ~isvector(z_profile) || length(r_profile) ~= length(z_profile)
        error('r_profile and z_profile must be vectors of the same length.');
    end
    if ~isscalar(thickness) || thickness <= 0
        error('thickness must be a positive scalar.');
    end
    if ~ischar(filename) && ~isstring(filename) || isempty(filename)
       error('filename must be a non-empty character vector or string.');
    end
    if any(r_profile < 0)
        error('r_profile coordinates must be non-negative.');
    end

    r_profile = r_profile(:); % Ensure column vector
    z_profile = z_profile(:); % Ensure column vector
    N_z = length(z_profile);

    if N_z < 2
        error('Profile must have at least 2 points.');
    end

    % --- Generate Inner Profile ---
    r_inner_profile = r_profile - thickness;

    % Check for negative or zero inner radius and clip if necessary
    if any(r_inner_profile <= 0)
        warning('Thickness is large relative to r_profile, resulting in zero or negative inner radius. Clipping inner radius at a small value (1e-6).');
        r_inner_profile = max(r_inner_profile, 1e-6); % Clip at small positive value
    end
    % Check if inner radius somehow became larger (shouldn't happen with positive thickness)
    if any(r_inner_profile >= r_profile)
         error('Calculated inner radius is not smaller than outer radius. Check profile points and thickness.');
    end

    % --- Generate Vertices ---
    % Create N_theta unique angles for the segments
    theta = linspace(0, 2*pi, n_revolution + 1); % N+1 points for N segments
    theta = theta(1:end-1); % Use N unique angles [0, 2pi)
    N_theta = length(theta); % Number of unique angular points = n_revolution

    % Calculate vertex coordinates using broadcasting
    % Outer surface points (Nz x N_theta matrices)
    X_outer = r_profile * cos(theta);
    Y_outer = r_profile * sin(theta);
    Z_outer = z_profile * ones(1, N_theta);

    % Inner surface points (Nz x N_theta matrices)
    X_inner = r_inner_profile * cos(theta);
    Y_inner = r_inner_profile * sin(theta);
    Z_inner = z_profile * ones(1, N_theta); % Z is the same for inner/outer

    % Combine vertices into one array [x, y, z]
    % Order: All outer points (z varies fastest, then theta), then all inner points
    V_outer = [X_outer(:), Y_outer(:), Z_outer(:)];
    V_inner = [X_inner(:), Y_inner(:), Z_inner(:)];
    V = [V_outer; V_inner];

    num_outer_verts_total = size(V_outer, 1); % Should be N_z * N_theta

    % Helper function to get linear index for vertices
    % i: index along profile (1 to N_z)
    % j: index around revolution (1 to N_theta)
    outer_idx = @(i, j) (j-1)*N_z + i;
    inner_idx = @(i, j) num_outer_verts_total + (j-1)*N_z + i;

    % --- Generate Faces (Triangles) ---
    num_faces_estimate = (N_z-1)*N_theta*2 * 2 + N_theta*2 * 2; % Outer wall + Inner wall + Top cap + Bottom cap
    F = zeros(num_faces_estimate, 3); % Preallocate face array
    f_count = 0; % Face counter

    for j = 1:N_theta
        j_next = mod(j, N_theta) + 1; % Next angular index, wraps around

        % 1. Outer Surface Wall (N_z-1 strips along Z)
        for i = 1:N_z-1
            % Define vertices of the quad for this segment
            p1 = outer_idx(i,   j);      % bottom-left
            p2 = outer_idx(i+1, j);      % top-left
            p3 = outer_idx(i+1, j_next); % top-right
            p4 = outer_idx(i,   j_next); % bottom-right

            % Create two triangles (ensure CCW order when viewed from outside)
            f_count = f_count + 1; F(f_count, :) = [p1, p2, p3];
            f_count = f_count + 1; F(f_count, :) = [p1, p3, p4];
        end

        % 2. Inner Surface Wall (N_z-1 strips along Z)
        for i = 1:N_z-1
            % Define vertices of the quad for this segment
            p1 = inner_idx(i,   j);      % bottom-left
            p2 = inner_idx(i+1, j);      % top-left
            p3 = inner_idx(i+1, j_next); % top-right
            p4 = inner_idx(i,   j_next); % bottom-right

            % Create two triangles (REVERSE order for CCW when viewed from outside -> normal points inwards)
            f_count = f_count + 1; F(f_count, :) = [p1, p3, p2]; % Reversed
            f_count = f_count + 1; F(f_count, :) = [p1, p4, p3]; % Reversed
        end

        % 3. Bottom Annulus (End Cap at i=1)
        % Define vertices of the annular segment
        p1_out = outer_idx(1, j);
        p2_in  = inner_idx(1, j);
        p3_in  = inner_idx(1, j_next);
        p4_out = outer_idx(1, j_next);

        % Create two triangles (ensure CCW order for normal pointing -Z)
        f_count = f_count + 1; F(f_count, :) = [p1_out, p3_in, p2_in]; % Reversed order for -Z normal
        f_count = f_count + 1; F(f_count, :) = [p1_out, p4_out, p3_in]; % Reversed order for -Z normal

        % 4. Top Annulus (End Cap at i=N_z)
        % Define vertices of the annular segment
        p1_out = outer_idx(N_z, j);
        p2_in  = inner_idx(N_z, j);
        p3_in  = inner_idx(N_z, j_next);
        p4_out = outer_idx(N_z, j_next);

        % Create two triangles (ensure CCW order for normal pointing +Z)
        f_count = f_count + 1; F(f_count, :) = [p1_out, p2_in, p3_in]; % Standard order for +Z normal
        f_count = f_count + 1; F(f_count, :) = [p1_out, p3_in, p4_out]; % Standard order for +Z normal
    end

    % Trim face array if preallocation was too large (shouldn't happen if logic is correct)
    if f_count < num_faces_estimate
        F = F(1:f_count, :);
    elseif f_count > num_faces_estimate
         warning('Actual face count (%d) exceeded estimate (%d). Check logic.', f_count, num_faces_estimate);
         % This indicates an error in the estimation or generation loop
    end


    % --- Write STL File ---
    try
        % Use the stlwrite function. Syntax might vary slightly.
        % Common syntaxes:
        % 1) stlwrite(filename, struct('faces', F, 'vertices', V))
        % 2) stlwrite(filename, F, V)
        % 3) stlwrite(filename, 'Faces', F, 'Vertices', V) % Newer syntax

        % Attempting newer Name-Value pair syntax first:
        stlwrite(filename, 'Faces', F, 'Vertices', V, 'Mode', 'binary'); % Binary is usually preferred
        fprintf('Successfully wrote binary STL file: %s\n', filename);

    catch ME
        % Handle cases where stlwrite might not exist or syntax differs
        if strcmp(ME.identifier, 'MATLAB:UndefinedFunction') || contains(ME.message, 'stlwrite')
             error(['stlwrite error: ', ME.message]);
        else
            % Rethrow other unexpected errors
            rethrow(ME);
        end
    end

end % function generate_revolution_stl