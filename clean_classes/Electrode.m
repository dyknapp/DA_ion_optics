classdef Electrode
    % Electrode class that builds an electrode profile from an external
    % function. The external function should return [rs, zs] given a
    % resolution (and optionally extra parameters).
    
    properties
        z_offset    % Position of FRONT EDGE along optical axis.
        rs          % 1xN array of radial coordinates
        zs          % 1xN array of axial coordinates
        segments    % (N-1)x4 array; each row is [r_start, z_start, r_end, z_end]
        profileFcn  % Function handle for generating the profile
        res         % Resolution overall used in the profile
        extraParams % Optional cell array of extra parameters for the profile function
    end
    
    methods
        function obj = Electrode(z_offset, profileFcn, res, varargin)
            % Constructor accepts a profile function and resolution. Extra
            % parameters for the profile function can be passed via varargin.
            %
            % Example usage:
            %   elec = Electrode(@ashfold_repeller_profile, 1);
            %
            %   or with extra parameters:
            %   elec = Electrode(@some_profile, 1, param1, param2);
            
            if ~isa(profileFcn, 'function_handle')
                error('profileFcn must be a function handle.');
            end
            obj.profileFcn = profileFcn;
            obj.res = res;
            obj.extraParams = varargin;
            
            % Generate the boundary points using the provided profile function.
            [obj.rs, obj.zs] = obj.profileFcn(obj.res, obj.extraParams{:});
            obj.zs = obj.zs + z_offset;
            
            % Compute the line segments from the boundary points.
            obj = obj.computeSegments();
        end
        
        function obj = computeSegments(obj)
            % Computes the line segments from consecutive (rs, zs) points.
            n = numel(obj.rs);
            segs = zeros(n-1, 4);
            for i = 1:(n-1)
                segs(i,:) = [obj.rs(i), obj.zs(i), obj.rs(i+1), obj.zs(i+1)];
            end
            obj.segments = segs;
        end
        
        function plotBoundary(obj)
            plot(obj.zs, obj.rs, '.-k');
        end
    end
end
