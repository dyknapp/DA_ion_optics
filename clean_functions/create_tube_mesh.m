function [V_tube, F_tube] = create_tube_mesh(path_points, radius, n_segments)
%CREATE_TUBE_MESH Generates vertices and faces for a tube around a 3D path.
%   Uses Parallel Transport to create a smoother frame along the curve.
%
%   INPUTS:
%       path_points (N x 3 matrix): Cartesian coordinates of points along the path.
%       radius (scalar): Radius of the tube.
%       n_segments (scalar): Number of segments around the tube circumference.
%
%   OUTPUTS:
%       V_tube (M x 3 matrix): Vertices of the tube mesh.
%       F_tube (P x 3 matrix): Faces (triangles) of the tube mesh.

    N_path = size(path_points, 1);
    if N_path < 2
        error('Path must have at least 2 points.');
    end
    if radius <= 0
        error('Radius must be positive.');
    end
    n_segments = max(round(n_segments), 3); % Need at least 3 segments

    V_tube = zeros(N_path * n_segments, 3); % Preallocate vertices
    F_tube = zeros((N_path - 1) * n_segments * 2, 3); % Preallocate faces

    % --- Calculate Tangents ---
    Tangents = zeros(N_path, 3);
    Tangents(1, :) = path_points(2, :) - path_points(1, :); % Forward difference for first point
    for i = 2:N_path-1
        Tangents(i, :) = path_points(i+1, :) - path_points(i-1, :); % Central difference
    end
    Tangents(N_path, :) = path_points(N_path, :) - path_points(N_path-1, :); % Backward difference for last point

    % Normalize tangents
    for i = 1:N_path
        normT = norm(Tangents(i, :));
        if normT > eps
            Tangents(i, :) = Tangents(i, :) / normT;
        else
            % Handle zero-length segment tangent (e.g., duplicate points)
            if i > 1
                Tangents(i, :) = Tangents(i-1, :); % Use previous tangent
            else % This should not happen if N_path >= 2
                Tangents(i, :) = [1, 0, 0]; % Default if first segment is zero length
            end
        end
    end

    % --- Calculate Frame (Normal N, Binormal B) using Parallel Transport ---
    Normals = zeros(N_path, 3);
    Binormals = zeros(N_path, 3);

    % Initialize first frame
    T1 = Tangents(1, :);
    % Find an initial normal vector non-parallel to T1
    arb_vec = [0, 1, 0];
    if abs(dot(T1, arb_vec)) > 0.99
        arb_vec = [1, 0, 0];
    end
    N1 = cross(T1, arb_vec);
    N1 = N1 / (norm(N1) + eps);
    B1 = cross(T1, N1); % Already normalized if T1, N1 orthonormal
    Normals(1, :) = N1;
    Binormals(1, :) = B1;

    % Transport frame along the curve
    for i = 2:N_path
        T_prev = Tangents(i-1, :);
        T_curr = Tangents(i, :);
        N_prev = Normals(i-1, :);
        B_prev = Binormals(i-1, :);

        % Rotation axis and angle to align T_prev with T_curr
        axis_rot = cross(T_prev, T_curr);
        sin_angle = norm(axis_rot);
        cos_angle = dot(T_prev, T_curr);

        if sin_angle > eps % Avoid division by zero if tangents are parallel
            axis_rot = axis_rot / sin_angle; % Normalize rotation axis
            angle_rot = atan2(sin_angle, cos_angle); % Angle of rotation

            % Rotate previous Normal and Binormal using Rodrigues' rotation formula
            % R(v) = v*cos(a) + (k x v)*sin(a) + k*(k . v)*(1-cos(a))
            N_curr = N_prev * cos_angle + cross(axis_rot, N_prev) * sin_angle + axis_rot * dot(axis_rot, N_prev) * (1 - cos_angle);
            B_curr = B_prev * cos_angle + cross(axis_rot, B_prev) * sin_angle + axis_rot * dot(axis_rot, B_prev) * (1 - cos_angle);

        else % Tangents are parallel (or anti-parallel)
            N_curr = N_prev; % Keep the frame the same
            B_curr = B_prev;
            % If anti-parallel (cos_angle ~ -1), need to flip normal? Check required.
            if cos_angle < -0.99
               N_curr = -N_curr; % Flip normal if tangent reversed direction
               B_curr = -B_curr;
            end
        end
        Normals(i, :) = N_curr / (norm(N_curr)+eps); % Re-normalize for safety
        Binormals(i, :) = B_curr / (norm(B_curr)+eps); % Re-normalize for safety

         % Ensure orthogonality (optional, but good practice)
         Binormals(i,:) = cross(Tangents(i,:), Normals(i,:));
         Binormals(i,:) = Binormals(i,:) / (norm(Binormals(i,:))+eps);
         Normals(i,:) = cross(Binormals(i,:), Tangents(i,:)); % Recalculate N from T, B
         Normals(i,:) = Normals(i,:) / (norm(Normals(i,:))+eps);

    end

    % --- Generate Vertices and Faces ---
    angles = linspace(0, 2*pi, n_segments + 1);
    angles = angles(1:end-1); % Use n_segments unique angles

    for i = 1:N_path
        P_i = path_points(i, :);
        N_i = Normals(i, :);
        B_i = Binormals(i, :);

        % --- Create Circle Vertices ---
        base_idx = (i - 1) * n_segments;
        for j = 1:n_segments
            V_tube(base_idx + j, :) = P_i + radius * (N_i * cos(angles(j)) + B_i * sin(angles(j)));
        end

        % --- Create Faces (Connecting circle i-1 to circle i) ---
        if i > 1
            base_idx_prev = (i - 2) * n_segments;
            face_base_idx = (i - 2) * n_segments * 2;

            for j = 1:n_segments
                j_next = mod(j, n_segments) + 1; % Wrap around index

                p1 = base_idx_prev + j;
                p2 = base_idx + j;
                p3 = base_idx + j_next;
                p4 = base_idx_prev + j_next;

                % Ensure correct winding order for outward normals
                F_tube(face_base_idx + (j-1)*2 + 1, :) = [p1, p2, p3];
                F_tube(face_base_idx + (j-1)*2 + 2, :) = [p1, p3, p4];
            end
        end
    end % End loop over path points

end % function create_tube_mesh