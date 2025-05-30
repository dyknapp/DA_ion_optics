function generate_lens_trajectory_stl(stlParams, electrodes, z_trajectory, r_trajectory_all)
%generate_lens_trajectory_stl Generates STL file for electrodes and optionally a trajectory solid.
%   Outputs revolution geometry with an optional angular cutaway for viewing.
%   The trajectory is represented as a solid of revolution with r=0 as the inner boundary.
%
%   INPUTS:
%       stlParams: Struct with STL parameters (R_outer, num_angular_steps,
%                  cutaway_angle_degrees, output_folder, export_trajectory)
%       electrodes: Vector of Electrode objects defining the lens geometry.
%       z_trajectory: Vector of z-coordinates for trajectory points (e.g., z_boundaries).
%       r_trajectory_all: Matrix of r-coordinates (num_z x num_rays). The last column
%                         is used if export_trajectory is true.

    % --- Extract Parameters ---
    R_outer = stlParams.R_outer;
    num_angular_steps = stlParams.num_angular_steps;
    % --- Use Cutaway Angle ---
    if isfield(stlParams, 'cutaway_angle_degrees') && isnumeric(stlParams.cutaway_angle_degrees) && isscalar(stlParams.cutaway_angle_degrees)
        cutaway_angle_degrees = stlParams.cutaway_angle_degrees;
        if cutaway_angle_degrees < 0 || cutaway_angle_degrees >= 360
             warning('generate_lens_trajectory_stl:invalidCutaway', ...
                     'Cutaway angle (%.2f deg) must be in [0, 360). Using 0 (full revolution).', cutaway_angle_degrees);
             cutaway_angle_degrees = 0;
        end
    else
        warning('generate_lens_trajectory_stl:missingCutaway', ...
                'Cutaway angle not specified or invalid in stlParams. Using 0 (full revolution).');
        cutaway_angle_degrees = 0; % Default to full revolution
    end
    cutaway_radians = deg2rad(cutaway_angle_degrees);
    output_folder = stlParams.output_folder;
    export_trajectory = stlParams.export_trajectory;

    % --- Determine Filename based on Cutaway ---
    is_full_revolution = (cutaway_radians <= eps(2*pi)); % Check if effectively zero
    if is_full_revolution
        combined_stl_filename = fullfile(output_folder, 'combined_lens_trajectory_full.stl');
    else
        combined_stl_filename = fullfile(output_folder, sprintf('combined_lens_trajectory_cutaway%.0f.stl', cutaway_angle_degrees));
    end

    % --- Input Validation ---
    if ~exist('stlwrite', 'file')
        error('generate_lens_trajectory_stl:stlwriteNotFound', ...
              'The stlwrite function is required. Please ensure it is on the MATLAB path.');
    end
    % Check trajectory dimensions only if exporting
    if export_trajectory
        if isempty(r_trajectory_all)
            warning('generate_lens_trajectory_stl:emptyTrajectory', ...
                    'r_trajectory_all is empty. Skipping trajectory export.');
            export_trajectory = false;
        elseif (size(r_trajectory_all, 1) ~= numel(z_trajectory))
             warning('generate_lens_trajectory_stl:trajectorySizeMismatch', ...
                     'Trajectory z (%d) and r (%d) dimensions mismatch. Skipping trajectory export.', ...
                     numel(z_trajectory), size(r_trajectory_all, 1));
             export_trajectory = false;
        end
    end

    % --- Setup for STL Generation (with Cutaway) ---
    if is_full_revolution
        % Generate theta for a full [0, 2*pi) range with num_angular_steps segments
        fprintf('  Generating full revolution (Cutaway Angle = 0 deg).\n');
        theta = linspace(0, 2*pi, num_angular_steps + 1);
        theta = theta(1:end-1); % Use N unique angles [0, 2pi)
        num_angular_points = num_angular_steps; % N distinct points/angles
        num_angular_segments = num_angular_steps; % N segments connecting N points cyclically
    else
        % Generate theta for the partial revolution [0, 2*pi - cutaway]
        end_angle = 2*pi - cutaway_radians;
        fprintf('  Generating partial revolution with cutaway (Angle = %.1f deg, Range = [0, %.3f rad]).\n', cutaway_angle_degrees, end_angle);
        % We need N+1 points to define N segments for the curved part
        theta = linspace(0, end_angle, num_angular_steps + 1);
        num_angular_points = num_angular_steps + 1; % N+1 distinct points/angles
        num_angular_segments = num_angular_steps;   % N segments connecting N+1 points linearly
    end
    theta = theta(:); % Ensure column vector

    % Ensure output folder exists
    if ~exist(output_folder, 'dir')
       mkdir(output_folder);
       fprintf('  Created output directory: %s\n', output_folder);
    end

    % --- Initialization for Combined Geometry ---
    all_vertices = []; % Stores vertices from all components
    all_faces = [];    % Stores faces from all components (with offset indices)
    vertex_offset = 0; % Running count of vertices added so far

    % ==========================================================
    % Part 1: Loop Through Electrodes to Generate Blocks
    % ==========================================================
    num_electrodes = numel(electrodes);
    fprintf('  Processing %d electrode blocks for combined STL...\n', num_electrodes);
    for i = 1:num_electrodes
        fprintf('    Processing electrode block %d...\n', i);

        % --- Get Profile Points for Electrode i ---
        profile_z = electrodes(i).zs; % Axial coordinates of inner boundary
        r_inner = electrodes(i).rs;   % Radial coordinates of inner boundary
        profile_z = profile_z(:);   % Ensure column
        r_inner = r_inner(:);       % Ensure column

        % Basic validation for electrode profile
        if numel(profile_z) ~= numel(r_inner)
            warning('generate_lens_trajectory_stl:profileMismatch', ...
                    'Skipping electrode block %d: Mismatch between number of z (%d) and r (%d) points.', ...
                    i, numel(profile_z), numel(r_inner));
            continue;
        end
        num_profile_points = numel(profile_z);
        if num_profile_points < 2
            warning('generate_lens_trajectory_stl:insufficientPoints', ...
                    'Skipping electrode block %d: Needs at least 2 profile points (found %d).', ...
                    i, num_profile_points);
            continue;
        end
        if any(r_inner < 0)
             warning('generate_lens_trajectory_stl:negativeRadius', ...
                    'Electrode block %d has negative inner radius. Check electrode definition.', i);
             r_inner = max(r_inner, 0); % Clamp at zero for safety
        end

        % Define outer radius profile (constant R_outer)
        r_outer = ones(num_profile_points, 1) * R_outer;
        if any(r_inner >= r_outer)
             warning('generate_lens_trajectory_stl:radiusConflict', ...
                    'Inner radius meets or exceeds outer radius (%.2f) for electrode %d. Check R_outer or electrode shape.', R_outer, i);
        end

        % ==========================================================
        % Part 2: Generate Vertices and Faces for the Current Electrode Block
        % ==========================================================
        fprintf('      Generating vertices and faces for electrode block %d...\n', i);

        % --- Calculate 3D Coordinates for Inner and Outer Surfaces ---
        X_inner = r_inner .* cos(theta'); % (num_profile x num_angular_points)
        Y_inner = r_inner .* sin(theta'); % (num_profile x num_angular_points)
        Z_inner = profile_z .* ones(1, num_angular_points); % (num_profile x num_angular_points)

        X_outer = r_outer .* cos(theta'); % (num_profile x num_angular_points)
        Y_outer = r_outer .* sin(theta'); % (num_profile x num_angular_points)
        Z_outer = profile_z .* ones(1, num_angular_points); % Same Z profile

        % --- Arrange Vertices ---
        inner_vertices = [X_inner(:), Y_inner(:), Z_inner(:)];
        outer_vertices = [X_outer(:), Y_outer(:), Z_outer(:)];
        electrode_vertices = [inner_vertices; outer_vertices];
        num_block_vertices = size(electrode_vertices, 1); % 2 * num_profile_points * num_angular_points
        num_surface_vertices = num_block_vertices / 2;

        % --- Generate Faces ---
        num_quads_curved = (num_profile_points - 1) * num_angular_segments; % Inner + Outer
        num_quads_caps = num_angular_segments * 2; % Top + Bottom
        num_quads_cut = 0;
        if ~is_full_revolution
            num_quads_cut = (num_profile_points - 1) * 2; % Start + End cut faces
        end
        total_faces = (num_quads_curved + num_quads_caps + num_quads_cut) * 2; % 2 triangles per quad

        electrode_faces = zeros(total_faces, 3);
        face_idx = 1;

        % Helper function for vertex indexing within the block
        inner_v_idx = @(j, k) (k-1)*num_profile_points + j;
        outer_v_idx = @(j, k) num_surface_vertices + (k-1)*num_profile_points + j;

        % --- Generate Curved and Cap Faces ---
        for k = 1:num_angular_segments % Iterate through the segments
            if is_full_revolution
                k_next = mod(k, num_angular_points) + 1;
            else
                k_next = k + 1;
            end

            % --- Generate Faces along the Profile (Inner, Outer) ---
            for j = 1:(num_profile_points - 1)
                j_next = j + 1;
                % Inner Surface Quad Indices & Triangles
                p1_in = inner_v_idx(j, k); p2_in = inner_v_idx(j_next, k); p3_in = inner_v_idx(j_next, k_next); p4_in = inner_v_idx(j, k_next);
                electrode_faces(face_idx,:) = [p1_in, p2_in, p3_in]; face_idx = face_idx + 1;
                electrode_faces(face_idx,:) = [p1_in, p3_in, p4_in]; face_idx = face_idx + 1;
                % Outer Surface Quad Indices & Triangles
                p1_out = outer_v_idx(j, k); p2_out = outer_v_idx(j_next, k); p3_out = outer_v_idx(j_next, k_next); p4_out = outer_v_idx(j, k_next);
                electrode_faces(face_idx,:) = [p1_out, p3_out, p2_out]; face_idx = face_idx + 1;
                electrode_faces(face_idx,:) = [p1_out, p4_out, p3_out]; face_idx = face_idx + 1;
            end

            % --- Bottom Cap Faces (j=1) ---
            p1_in = inner_v_idx(1, k); p4_in = inner_v_idx(1, k_next); p1_out = outer_v_idx(1, k); p4_out = outer_v_idx(1, k_next);
            electrode_faces(face_idx,:) = [p1_in, p1_out, p4_out]; face_idx = face_idx + 1;
            electrode_faces(face_idx,:) = [p1_in, p4_out, p4_in]; face_idx = face_idx + 1;

            % --- Top Cap Faces (j=num_profile_points) ---
            j_top = num_profile_points;
            p2_in = inner_v_idx(j_top, k); p3_in = inner_v_idx(j_top, k_next); p2_out = outer_v_idx(j_top, k); p3_out = outer_v_idx(j_top, k_next);
            electrode_faces(face_idx,:) = [p2_in, p3_in, p2_out]; face_idx = face_idx + 1;
            electrode_faces(face_idx,:) = [p3_in, p3_out, p2_out]; face_idx = face_idx + 1;
        end % End loop over angular segments (k)

        % --- Generate Cut Faces (if not full revolution) ---
        if ~is_full_revolution
            fprintf('      Generating cut faces for electrode block %d...\n', i);
            k_start = 1; k_end = num_angular_points;
            for j = 1:(num_profile_points - 1)
                j_next = j + 1;
                % Start Cut Face (k=1)
                p1 = inner_v_idx(j, k_start); p2 = outer_v_idx(j, k_start); p3 = outer_v_idx(j_next, k_start); p4 = inner_v_idx(j_next, k_start);
                electrode_faces(face_idx,:) = [p1, p2, p3]; face_idx = face_idx + 1;
                electrode_faces(face_idx,:) = [p1, p3, p4]; face_idx = face_idx + 1;
                % End Cut Face (k=k_end)
                p1 = inner_v_idx(j, k_end); p2 = outer_v_idx(j, k_end); p3 = outer_v_idx(j_next, k_end); p4 = inner_v_idx(j_next, k_end);
                electrode_faces(face_idx,:) = [p1, p3, p2]; face_idx = face_idx + 1; % Flipped order for outward normal
                electrode_faces(face_idx,:) = [p1, p4, p3]; face_idx = face_idx + 1; % Flipped order
            end
        end % End if ~is_full_revolution

         if face_idx - 1 ~= total_faces % Check face count
             warning('generate_lens_trajectory_stl:faceCountMismatch', ...
                     'Electrode Block %d: Expected %d faces, generated %d faces.', i, total_faces, face_idx - 1);
             electrode_faces = electrode_faces(1:face_idx-1, :); % Trim if necessary
         end

        % --- Append to Global Lists ---
        fprintf('      Appending %d vertices and %d faces for electrode block %d.\n', num_block_vertices, size(electrode_faces, 1), i);
        all_vertices = [all_vertices; electrode_vertices];
        electrode_faces_offset = electrode_faces + vertex_offset;
        all_faces = [all_faces; electrode_faces_offset];
        vertex_offset = vertex_offset + num_block_vertices;

    end % End loop over electrodes

    % ==========================================================
    % Part 2.5: Generate Trajectory Solid (if requested)
    % ==========================================================
    if export_trajectory
        fprintf('  Processing trajectory solid export...\n');

        % --- Get Trajectory Profile ---
        traj_profile_z = z_trajectory(:);
        % Use last ray, ensure non-negative radius
        traj_r_outer = max(r_trajectory_all(:, end), 0);
        traj_r_inner = zeros(size(traj_profile_z)); % Inner radius is zero (on axis)
        num_traj_profile_points = numel(traj_profile_z);

        if num_traj_profile_points < 2
            warning('generate_lens_trajectory_stl:insufficientTrajPoints', ...
                    'Skipping trajectory solid export: Needs at least 2 trajectory points (found %d).', ...
                     num_traj_profile_points);
        else
            fprintf('    Generating vertices and faces for trajectory solid...\n');

            % --- Calculate Trajectory 3D Coordinates ---
            % Inner coordinates will be on the Z-axis (X=0, Y=0)
            Traj_X_inner = traj_r_inner .* cos(theta'); % All zeros
            Traj_Y_inner = traj_r_inner .* sin(theta'); % All zeros
            Traj_Z_inner = traj_profile_z .* ones(1, num_angular_points);
            % Outer coordinates follow the trajectory radius
            Traj_X_outer = traj_r_outer .* cos(theta');
            Traj_Y_outer = traj_r_outer .* sin(theta');
            Traj_Z_outer = traj_profile_z .* ones(1, num_angular_points); % Same Z profile

            % --- Arrange Trajectory Vertices ---
            traj_inner_vertices = [Traj_X_inner(:), Traj_Y_inner(:), Traj_Z_inner(:)];
            traj_outer_vertices = [Traj_X_outer(:), Traj_Y_outer(:), Traj_Z_outer(:)];
            traj_vertices = [traj_inner_vertices; traj_outer_vertices];
            num_traj_block_vertices = size(traj_vertices, 1);
            num_traj_surface_vertices = num_traj_block_vertices / 2;

            % --- Generate Trajectory Faces ---
            % Faces needed: Outer surface, End Caps (simplified), Cut Faces (if applicable)
            num_quads_outer_traj = (num_traj_profile_points - 1) * num_angular_segments;
            num_tris_caps_traj = num_angular_segments * 2; % Simplified caps (1 triangle each end per segment)
            num_quads_cut_traj = 0;
            if ~is_full_revolution
                num_quads_cut_traj = (num_traj_profile_points - 1) * 2; % Start + End cut faces
            end
            total_traj_faces = num_quads_outer_traj * 2 + num_tris_caps_traj + num_quads_cut_traj * 2;

            traj_faces = zeros(total_traj_faces, 3);
            traj_face_idx = 1;

            % Helper function for trajectory vertex indexing
            traj_inner_v_idx = @(j, k) (k-1)*num_traj_profile_points + j;
            traj_outer_v_idx = @(j, k) num_traj_surface_vertices + (k-1)*num_traj_profile_points + j;

            % --- Generate Outer Surface and Simplified Cap Faces ---
            for k = 1:num_angular_segments % Iterate through the segments
                if is_full_revolution
                    k_next = mod(k, num_angular_points) + 1;
                else
                    k_next = k + 1;
                end

                % --- Generate Outer Surface Faces ---
                for j = 1:(num_traj_profile_points - 1)
                    j_next = j + 1;
                    p1_out = traj_outer_v_idx(j, k); p2_out = traj_outer_v_idx(j_next, k); p3_out = traj_outer_v_idx(j_next, k_next); p4_out = traj_outer_v_idx(j, k_next);
                    % Outer Triangles (Normals point outwards)
                    traj_faces(traj_face_idx,:) = [p1_out, p3_out, p2_out]; traj_face_idx = traj_face_idx + 1;
                    traj_faces(traj_face_idx,:) = [p1_out, p4_out, p3_out]; traj_face_idx = traj_face_idx + 1;
                end

                % --- Simplified Bottom Cap Faces (j=1) ---
                % Connect axis point to outer edge segment
                p_axis_k = traj_inner_v_idx(1, k); % Inner point at k (same for k_next)
                p_out_k = traj_outer_v_idx(1, k);
                p_out_k_next = traj_outer_v_idx(1, k_next);
                % Single triangle per segment (Normal approx -Z)
                traj_faces(traj_face_idx,:) = [p_axis_k, p_out_k, p_out_k_next]; traj_face_idx = traj_face_idx + 1;

                % --- Simplified Top Cap Faces (j=num_traj_profile_points) ---
                j_top = num_traj_profile_points;
                p_axis_k_top = traj_inner_v_idx(j_top, k); % Inner point at k (same for k_next)
                p_out_k_top = traj_outer_v_idx(j_top, k);
                p_out_k_next_top = traj_outer_v_idx(j_top, k_next);
                 % Single triangle per segment (Normal approx +Z)
                traj_faces(traj_face_idx,:) = [p_axis_k_top, p_out_k_next_top, p_out_k_top]; traj_face_idx = traj_face_idx + 1;

            end % End loop over angular segments (k)

            % --- Generate Trajectory Cut Faces (if not full revolution) ---
            if ~is_full_revolution
                fprintf('      Generating cut faces for trajectory solid...\n');
                k_start = 1; k_end = num_angular_points;
                for j = 1:(num_traj_profile_points - 1)
                    j_next = j + 1;
                    % Start Cut Face (k=1) - Connects axis to outer edge
                    p1 = traj_inner_v_idx(j, k_start); p2 = traj_outer_v_idx(j, k_start); p3 = traj_outer_v_idx(j_next, k_start); p4 = traj_inner_v_idx(j_next, k_start);
                    traj_faces(traj_face_idx,:) = [p1, p2, p3]; traj_face_idx = traj_face_idx + 1;
                    traj_faces(traj_face_idx,:) = [p1, p3, p4]; traj_face_idx = traj_face_idx + 1;
                    % End Cut Face (k=k_end) - Connects axis to outer edge
                    p1 = traj_inner_v_idx(j, k_end); p2 = traj_outer_v_idx(j, k_end); p3 = traj_outer_v_idx(j_next, k_end); p4 = traj_inner_v_idx(j_next, k_end);
                    traj_faces(traj_face_idx,:) = [p1, p3, p2]; traj_face_idx = traj_face_idx + 1; % Flipped order
                    traj_faces(traj_face_idx,:) = [p1, p4, p3]; traj_face_idx = traj_face_idx + 1; % Flipped order
                end
            end % End if ~is_full_revolution

            if traj_face_idx - 1 ~= total_traj_faces % Check face count
                warning('generate_lens_trajectory_stl:trajFaceCountMismatch', ...
                        'Trajectory Solid: Expected %d faces, generated %d faces.', total_traj_faces, traj_face_idx - 1);
                traj_faces = traj_faces(1:traj_face_idx-1, :); % Trim if necessary
            end

            % --- Append Trajectory Geometry to Global Lists ---
            fprintf('      Appending %d vertices and %d faces for trajectory solid.\n', num_traj_block_vertices, size(traj_faces, 1));
            all_vertices = [all_vertices; traj_vertices];
            traj_faces_offset = traj_faces + vertex_offset; % Apply current offset
            all_faces = [all_faces; traj_faces_offset];
            vertex_offset = vertex_offset + num_traj_block_vertices; % Update offset

        end % End if num_traj_profile_points >= 2
    end % End if export_trajectory

    % ==========================================================
    % Part 3: Final Export of Combined Geometry
    % ==========================================================
    fprintf('  Exporting combined geometry...\n');
    fprintf('    Total combined vertices: %d\n', size(all_vertices, 1));
    fprintf('    Total combined faces: %d\n', size(all_faces, 1));

    if isempty(all_vertices) || isempty(all_faces)
        warning('generate_lens_trajectory_stl:noGeometry', ...
                'No geometry generated (empty vertices or faces). Skipping final STL export.');
    else
        % --- Create Combined Triangulation Object ---
        fprintf('    Creating final triangulation object...\n');
        try
            T_combined = triangulation(all_faces, all_vertices);
            fprintf('    Triangulation object created successfully.\n');

            % --- Check for Unreferenced Vertices (Debugging) ---
            try
                referenced_verts_idx = unique(T_combined.ConnectivityList(:));
                if ~isempty(referenced_verts_idx) % Check if ConnectivityList is not empty
                    all_verts_idx = (1:size(all_vertices, 1))';
                    unreferenced_verts_idx = setdiff(all_verts_idx, referenced_verts_idx);
                    if ~isempty(unreferenced_verts_idx)
                        fprintf('    WARNING: Found %d unreferenced vertices.\n', numel(unreferenced_verts_idx));
                    else
                        fprintf('    No unreferenced vertices found.\n');
                    end
                else
                     fprintf('    WARNING: ConnectivityList is empty, cannot check for unreferenced vertices.\n');
                end
            catch ME_debug
                 warning('generate_lens_trajectory_stl:debugFail', ...
                         'Could not perform check for unreferenced vertices. Error: %s', ME_debug.message);
            end

            % --- Export Combined Geometry to STL ---
            fprintf('    Writing combined STL file: %s...\n', combined_stl_filename);
            try
                stlwrite(T_combined, combined_stl_filename, 'binary');
                fprintf('    Successfully wrote combined binary STL file.\n');
            catch ME_write
                 warning('generate_lens_trajectory_stl:stlWriteFail', ...
                         'Failed to write combined STL file %s. Error: %s', combined_stl_filename, ME_write.message);
                 % Fallback to ASCII
                 try
                    fprintf('    Attempting fallback to ASCII mode...\n');
                    stlwrite(T_combined, combined_stl_filename, 'mode', 'ascii');
                    fprintf('    Successfully wrote combined ASCII STL file as fallback.\n');
                 catch ME_write_ascii
                    warning('generate_lens_trajectory_stl:stlWriteAsciiFail', 'Failed to write ASCII STL as well: %s', ME_write_ascii.message);
                 end
            end % End try-catch stlwrite

        catch ME_tri
            warning('generate_lens_trajectory_stl:triangulationFail', ...
                    'Failed to create final combined triangulation object. Error: %s\n Check for issues like non-manifold geometry or inconsistent face definitions.', ...
                    ME_tri.message);
            % Fallback: try writing directly from faces/vertices
            fprintf('    Attempting direct STL write from vertices/faces...\n');
             try
                stlwrite(combined_stl_filename, 'Faces', all_faces, 'Vertices', all_vertices, 'mode', 'binary');
                fprintf('    Successfully wrote combined binary STL file directly.\n');
            catch ME_write_direct
                 warning('generate_lens_trajectory_stl:stlWriteDirectFail', ...
                         'Direct STL write also failed (binary). Error: %s', ME_write_direct.message);
                 % Fallback to ASCII
                 try
                    fprintf('    Attempting direct fallback to ASCII mode...\n');
                    stlwrite(combined_stl_filename, 'Faces', all_faces, 'Vertices', all_vertices, 'mode', 'ascii');
                    fprintf('    Successfully wrote combined ASCII STL file directly as fallback.\n');
                 catch ME_write_direct_ascii
                    warning('generate_lens_trajectory_stl:stlWriteDirectAsciiFail', 'Direct ASCII STL write also failed: %s', ME_write_direct_ascii.message);
                 end
            end
        end % End try-catch triangulation
    end % End if isempty check

    fprintf('  STL generation process finished.\n');

end % function generate_lens_trajectory_stl
