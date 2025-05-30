function [electrodes, z_boundaries_elements] = createElectrodeStack(electrodeData, opts)
    % createElectrodeStack Creates a vector of Electrode objects with computed z positions,
    % allowing for explicit spacing commands.
    %
    % INPUTS:
    %   electrodeData - A cell array where each element is either:
    %                   1) A structure defining an electrode:
    %                      .shape     - A function handle for the profile.
    %                      .thickness - (Optional) Scalar thickness (default=0).
    %                      .extraParams - (Optional) Cell array of extra parameters.
    %                   2) A cell array {'space', gap_value} to add explicit spacing.
    %
    %   opts          - A structure with optional fields:
    %                      .z0      - Starting z position (default = 0). Defines the position
    %                                 before the first spacing/element.
    %                      .spacing - Default gap added *before* each electrode struct (default = 1).
    %                      .res     - Overall resolution for electrode profiles (default = 1024).
    %
    % OUTPUT:
    %   electrodes    - A column vector of Electrode objects created.
    %   z_boundaries_elements - Vector of z-coordinates at element boundaries.

    arguments
        electrodeData {mustBeA(electrodeData, 'cell')} % Positional argument 1
        opts struct = struct() % Positional argument 2: opts struct (optional, defaults to empty)
    end

    if ~isfield(opts, 'z0') || isempty(opts.z0) || ~isscalar(opts.z0) || ~isnumeric(opts.z0)
        opts_z0 = 0; % Default z0
    else
        opts_z0 = opts.z0;
    end
    if ~isfield(opts, 'spacing') || isempty(opts.spacing) || ~isscalar(opts.spacing) || ~isnumeric(opts.spacing) || opts.spacing < 0
        opts_spacing = 1; % Default spacing
    else
        % Ensure non-negative if provided
        opts_spacing = max(opts.spacing, 0);
    end
    if ~isfield(opts, 'res') || isempty(opts.res) || ~isscalar(opts.res) || ~isnumeric(opts.res) || opts.res <= 0
        opts_res = 1024; % Default res
    else
        % Ensure positive if provided
         opts_res = max(opts.res, 1); % Or some small positive number
    end

    n_items = numel(electrodeData);
    electrode_objects_list = cell(n_items, 1); % Preallocate list for Electrode objects
    num_electrodes_created = 0;
    z_boundaries_elements = zeros(n_items + 1, 1); % Initialize boundary storage
    current_z = opts.z0; % Start at the initial z position
    z_boundaries_elements(1) = current_z; % Store initial boundary

    fprintf('Building electrode stack:\n');
    fprintf('  Initial z = %.3f\n', current_z);

    for i = 1:n_items
         item = electrodeData{i};

         if isstruct(item) % --- Handle Electrode Definition ---
             data = item;

             % 1. Add default spacing BEFORE the electrode
             current_z = current_z + opts.spacing;
             fprintf('  Item %d (Electrode): Adding default spacing %.3f -> z = %.3f\n', i, opts.spacing, current_z);

             % 2. Determine electrode thickness
             if isfield(data, 'thickness') && isscalar(data.thickness) && isnumeric(data.thickness) && data.thickness >= 0
                 thickness = data.thickness;
             else
                 thickness = 0; % Default thickness if not specified or invalid
                 if isfield(data, 'thickness')
                     warning('createElectrodeStack:InvalidThickness', 'Invalid thickness provided for item %d. Using thickness=0.', i);
                 end
             end

             % 3. Set z_offset (front edge position)
             z_offset = current_z;
             fprintf('    Creating Electrode at front edge z = %.3f (thickness=%.3f)\n', z_offset, thickness);


             % 4. Validate shape function handle
             if ~isfield(data, 'shape') || ~isa(data.shape, 'function_handle')
                 warning('createElectrodeStack:InvalidShape', 'Skipping item %d: Missing or invalid shape function handle.', i);
                 % Advance z by thickness anyway to maintain position for next element?
                 % Let's advance z by thickness, assuming space was allocated.
                  current_z = current_z + thickness;
                  fprintf('    Skipped electrode creation, advancing z by thickness -> z = %.3f\n', current_z);
                 continue; % Skip to next item
             end

             % 5. Create the Electrode object
             try
                 if isfield(data, 'extraParams') && iscell(data.extraParams)
                     new_electrode = Electrode(z_offset, data.shape, opts.res, data.extraParams{:});
                 else
                     new_electrode = Electrode(z_offset, data.shape, opts.res);
                 end
                 % Store the created electrode
                 num_electrodes_created = num_electrodes_created + 1;
                 electrode_objects_list{num_electrodes_created} = new_electrode;
             catch ME
                  warning('createElectrodeStack:ElectrodeError', 'Error creating Electrode object for item %d: %s. Skipping item.', i, ME.message);
                  % Advance z by thickness as if space was allocated
                  current_z = current_z + thickness;
                  fprintf('    Skipped electrode creation due to error, advancing z by thickness -> z = %.3f\n', current_z);
                  continue; % Skip to next item
             end

             % 6. Update current_z to the position AFTER this electrode
             current_z = current_z + thickness;
             fprintf('    Electrode ends at z = %.3f\n', current_z);


         elseif iscell(item) && numel(item) == 2 && strcmpi(item{1}, 'space') % --- Handle Space Command ---
             gap_value = item{2};
             if isnumeric(gap_value) && isscalar(gap_value)
                 current_z = current_z + gap_value;
                 fprintf('  Item %d (Space): Adding explicit space %.3f -> z = %.3f\n', i, gap_value, current_z);
             else
                  warning('createElectrodeStack:InvalidSpace', 'Invalid gap value for item %d. Ignoring space command.', i);
             end

         else % --- Handle Invalid Item Type ---
             warning('createElectrodeStack:InvalidItem', 'Item %d is not a valid electrode struct or a {''space'', value} cell. Ignoring.', i);
         end
             z_boundaries_elements(i+1) = current_z; % Store end position for this item
    end % End loop through electrodeData

    % Remove empty cells from the list
    electrode_objects_list = electrode_objects_list(1:num_electrodes_created);

    if num_electrodes_created > 0
        % Convert the cell array into a column vector of Electrode objects.
        electrodes = vertcat(electrode_objects_list{:});
        fprintf('Successfully created %d Electrode objects.\n', num_electrodes_created);
    else
        electrodes = Electrode.empty(0,1); % Return empty Electrode array
        fprintf('No Electrode objects were created.\n');
    end

end

% Helper function for argument validation (can be kept or removed if not needed elsewhere)
function mustBeA(x, typeName)
    if ~isa(x, typeName) % Use standard isa function
        error('Input must be of type %s.', typeName);
    end
end
function mustBeNonnegative(x)
    if ~isnumeric(x) || ~isscalar(x) || x < 0
        error('Value must be a non-negative scalar number.');
    end
end
function mustBePositive(x)
     if ~isnumeric(x) || ~isscalar(x) || x <= 0
        error('Value must be a positive scalar number.');
    end
end