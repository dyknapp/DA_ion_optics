function x_opt = optimize_einzel_lens()
%OPTIMIZE_EINZEL_LENS  Optimize the six lens voltages to minimize RMS spot size.
%
%   x_opt = OPTIMIZE_EINZEL_LENS() varies
%       [el1, el2, dl1, dl2, dl3, dl4] ∈ [–5, 0.5]
%   to minimize the RMS radius of the exiting ray‐fan.  You can hit
%   Ctrl+C at any time; the best‐so‐far solution will be returned.

  %% ----------------------
  % 1) GLOBAL SIMULATION SETUP
  %% ----------------------
  simParams.res       = 1;
  simParams.spacing   = 0;
  simParams.z0        = 0;
  simParams.R         = 1;

  daParams.daOrder    = 5;
  daParams.nVars      = 2;
  daParams.q_charge   = 1.0;
  daParams.m_mass     = 1.0;
  daParams.E0         = 0.6;

  dopriParams.z_start_overall = 30;
  dopriParams.z_end_overall   = 462;

  integratorOpts = struct( ...
    'RelTol',    1e-8, ...
    'AbsTol',    1e-11, ...
    'Stats',     'off', ...
    'MaxStep',   2.0);

% --- Electrode Geometry Definition (Uses createElectrodeStack format) ---
% Define electrode shapes and thicknesses
BENDER_PLACEHOLDER = 0.5 * 43 * pi + 7.95*2;

% Deccelerator stack params.
tt = 2.0;
pt = 4.0;
ids = {15.5, 20, 17, 16, 6, 32};
ls  = {33.0, 46, 12, 12, 6, 26};
ts  = {33 + 3, 46 + 2, 12 + 2,  12 + 2, 6 + 2, 26 + 3};

% Actual electrode data
electrodeData = {
    struct('thickness', 12, 'shape', @(res) ETH_endcap(res))            % 1
    struct('thickness', 12, 'shape', @(res) ETH_profile(res))           % 2
    struct('thickness', 12, 'shape', @(res) ETH_profile(res))           % 3
    % INTERACTION REGION
    struct('thickness', 12, 'shape', @(res) ETH_profile(res))           % 4
    struct('thickness', 12, 'shape', @(res) ETH_profile(res))           % 5
    struct('thickness', 12, 'shape', @(res) ETH_narrow(res))            % 6
    struct('thickness', 12, 'shape', @(res) ETH_profile(res))           % 7
    struct('thickness', 24, 'shape', @(res) ashfold_steer_tube(res))    % 8
    struct('thickness', 12, 'shape', @(res) ETH_profile(res))           % 9
    struct('thickness', 12, 'shape', @(res) ashfold_empty_tube(res))    % 10
    struct('thickness', 12, 'shape', @(res) ETH_endcap_right(res))      % 11
    {'space', BENDER_PLACEHOLDER}                                       % 12: Bender Placeholder
    struct('thickness', 12, 'shape', @(res) ETH_endcap(res))            % 13
    struct('thickness', 24, 'shape', @(res) ashfold_steer_tube(res))    % 14
    struct('thickness', 12, 'shape', @(res) ETH_endcap_right(res))      % 15
    {'space', 37.3}
    struct('thickness', ts{1}, 'shape', ...
        @(res) tube_in_plate(ids{1}, tt, ls{1}, pt, 37.5, res))
    struct('thickness', ts{2}, 'shape', ...
        @(res) tube_in_plate(ids{2}, tt, ls{2}, pt, 37.5, res))
    struct('thickness', ts{3}, 'shape', ...
        @(res) tube_in_plate(ids{3}, tt, ls{3}, pt, 37.5, res))
    struct('thickness', ts{4}, 'shape', ...
        @(res) tube_in_plate(ids{4}, tt, ls{4}, pt, 37.5, res))
    struct('thickness', ts{5}, 'shape', ...
        @(res) tube_in_plate(ids{5}, tt, ls{5}, pt, 37.5, res))
    struct('thickness', ts{6}, 'shape', ...
        @(res) tube_in_plate(ids{6}, tt, ls{6}, pt, 37.5, res))
    struct('thickness', 8, 'shape', ...
        @(res) tube_in_plate(3, 0.5, 6, 2, 37.5, res))
};

  % load pre‐computed quadrupole bender DA maps
  ld = load('MAT007_clean_DA/bender_segmented_data.mat', ...
            'transfer_maps','ref_states','z_coords');
  transfer_maps_bender = ld.transfer_maps;
  ref_states_bender    = ld.ref_states;
  z_coords_bender      = ld.z_coords;

  %% ----------------------
  % 3) BUILD ELECTRODE STACK & RUN BEM ONCE
  %% ----------------------
  opts_stack = struct('z0',simParams.z0,'spacing',simParams.spacing,'res',simParams.res);
  [electrodes_full_stack, z_boundaries] = ...
      createElectrodeStack(electrodeData, opts_stack);

  % locate bender insertion region
  bidx = find(cellfun(@(c)iscell(c)&&strcmpi(c{1},'space'), electrodeData),1);
  z_bend_start = z_boundaries(bidx) + (7.95 - (116 - 110.696));
  z_bend_end   = z_boundaries(bidx+1);
  z_coords_bender = z_coords_bender + z_bend_start;

  % one‐time BEM solve → qs & geometry table
  [qs, bemTable] = BEM_monopole(electrodes_full_stack, 1.0, false);

  % package static voltage‐independent params for DA integration
  bemParams.q_charge = daParams.q_charge;
  bemParams.m_mass   = daParams.m_mass;
  bemParams.E0       = daParams.E0;

  %% ----------------------
  % 4) PRECOMPUTE BENDER DA MAP SEQUENCE
  %% ----------------------
  sim_daOrder = daParams.daOrder;
  sim_nVars   = daParams.nVars;
  r_da  = DiffAlg.var(1,sim_daOrder,sim_nVars);
  pr_da = DiffAlg.var(2,sim_daOrder,sim_nVars);
  one_da= DiffAlg.one(sim_nVars,sim_daOrder);

  nseg = numel(transfer_maps_bender);
  S_maps_bender_segmented = cell(1,nseg);
  for i = 1:nseg
    M   = transfer_maps_bender{i};
    X0  = ref_states_bender(:,i);
    X1  = ref_states_bender(:,i+1);
    Off = X1 - M*X0;
    S_maps_bender_segmented{i} = { ...
      r_da*M(1,1) + pr_da*M(1,2) + one_da*Off(1), ...
      r_da*M(2,1) + pr_da*M(2,2) + one_da*Off(2)  ...
    };
  end

  %% ----------------------
  % 5) RAY‐TRACE & OPTIMIZER SETTINGS
  %% ----------------------
  traceParams.num_r_rays = 16;
  traceParams.num_s_rays = 16;
  traceParams.max_r0_dev = 0.5;
  traceParams.max_s0_dev = 0.001;
  traceParams.center_r0  = 0.0;
  traceParams.center_s0  = 0.0;

  visParams.electrodes_full_stack = electrodes_full_stack;
  visParams.simParams             = simParams;
  trace_eval_order = sim_daOrder;

  %% ----------------------
  % 6) RUN FMINCON WITH INTERRUPT HANDLING
  %% ----------------------
  x0 = [0.1; 0.2; 0.3; 0.4; 0.3; 0.2; -0.1; -0.2;];
  LB = -2*ones(8,1);
  UB =  0.5*ones(8,1);

  opts = optimoptions('fmincon', ...
           'Display','iter', ...
           'MaxFunctionEvaluations', 5000, ...
           'UseParallel', false, ...
           'Algorithm', 'sqp', ...
           'OutputFcn', @outfun);

  % container for best‐so‐far
  bestSolution.x    = x0;
  bestSolution.fval = Inf;

  try
    [x_opt, fval] = fmincon(@objective, x0, [],[],[],[], LB, UB, [], opts);
  catch ME
    if strcmp(ME.identifier,'MATLAB:execution:Interrupted')
      fprintf('Optimization interrupted—returning best‐so‐far (RMS=%.6g).\n', ...
              bestSolution.fval);
      x_opt = bestSolution.x;
      fval  = bestSolution.fval;
    else
      rethrow(ME);
    end
  end

  fprintf('\nFinal RMS spot size = %.6f\n', fval);
  fprintf('Optimal voltages: el1=%.4f, el2=%.4f, dl1=%.4f, dl2=%.4f, dl3=%.4f, dl4=%.4f, dl5=%.4f, ecv=%.4f\n', ...
          x_opt);

  if nargout<1, clear x_opt; end

  %% ----------------------
  % Nested OUTPUTFCN: save best‐so‐far on each iter
  %% ----------------------
    function stop = outfun(x,optimValues,state)
    stop = false;
    switch state
      case 'init'
        bestSolution.x    = x;
        bestSolution.fval = optimValues.fval;
        save('best_solution.mat','bestSolution');
      case 'iter'
        if optimValues.fval < bestSolution.fval
          bestSolution.x    = x;
          bestSolution.fval = optimValues.fval;
          save('best_solution.mat','bestSolution');
        end
      case 'done'
        % nothing
    end
  end

  %% ----------------------
  % Nested OBJECTIVE: update voltages, integrate, trace, compute RMS
  %% ----------------------
  function rmsSize = objective(x)
    el1 = x(1);  se1 = x(2);
    se2 = x(3);  dl1 = x(4);
    dl2 = x(5);  dl3 = x(6);
    dl4 = x(7);  ecv = x(8);

    % assemble Vs and update bemParams
    %
    Vs = [1.0; 0.75; 0.50; 0.25; ...    % 1-4: uniform field tube
          0.0; el1; 0.0; ...            % 5-7: Einzel 1
          se1; 0.0; 0.0; 0.0; ...       % 8-11: steering electrode 1
          % No voltage for the bender 'space' (item 10)
          0.0; se2; 0.0;...             % 10-12: Einzel 2
          0.0; dl1; dl2; dl3; dl4; ecv; % Decelerator assembly
          ecv; % endcap electrode
         ];
    bemParams.Vs  = Vs;
    bemParams.qVs = qs' * Vs;

    % integrate up to bender
    p1 = struct('z_start', dopriParams.z_start_overall, ...
                'z_end',   z_bend_start);
    [z1,S1,~] = run_dopri87_da_integration(p1, daParams, bemParams, bemTable, integratorOpts);

    % integrate after bender
    p2 = struct('z_start', z_bend_end, ...
                'z_end',   dopriParams.z_end_overall);
    [z2,S2,~] = run_dopri87_da_integration(p2, daParams, bemParams, bemTable, integratorOpts);

    % multi‐segment ray‐fan trace
    [rpos,~,~] = trace_ray_fan_multisegment_da_map( ...
                    traceParams, ...
                    struct('z_steps', z1, 'S_steps',{S1}), ...
                    S_maps_bender_segmented, ...
                    z_coords_bender, ...
                    struct('z_steps', z2, 'S_steps',{S2}), ...
                    visParams, ...
                    trace_eval_order, ...
                    true);
    drawnow;

    % extract final radii & compute RMS spot size
    rf = squeeze(rpos(end, :));
    rmsSize = sqrt(mean(rf(:).^2));
  end
end
