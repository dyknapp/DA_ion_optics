function x_opt = optimize_einzel_lens_sa()
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

  dopriParams.z_start_overall = 28.75;
  dopriParams.z_end_overall   = 477.7;

  integratorOpts = struct( ...
    'RelTol',    1e-8, ...
    'AbsTol',    1e-11, ...
    'Stats',     'off', ...
    'MaxStep',   2.0, ...
    'verbose',   false);


% Electrode geometry definition
BENDER_PLACEHOLDER = 0.5 * 43 * pi + 7.95*2;
tt = 2.0; pt = 4.0;
ids = {15.5,20,17,16,6,32};
ls  = {33.0,46,12,12,6,26};
ts  = {33,46+4.5,12+4.5,12+4.5,6+4.5,26};

electrodeData = { ...
struct('thickness',11.5,'shape',@(res)ETH_endcap(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_profile(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_profile(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_profile(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_profile(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_narrow(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_profile(res)), ...
struct('thickness',23.5,'shape',@(res)ashfold_steer_tube(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_profile(res)), ...
struct('thickness',11.5,'shape',@(res)ashfold_empty_tube(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_endcap_right(res)), ...
{'space', BENDER_PLACEHOLDER}, ...
struct('thickness',11.5,'shape',@(res)ETH_endcap(res)), ...
struct('thickness',23.5,'shape',@(res)ashfold_steer_tube(res)), ...
struct('thickness',11.5,'shape',@(res)ETH_endcap_right(res)), ...
struct('thickness',55.0,'shape',@(res) tube(25, 50, 32*res)), ...
struct('thickness',ts{1},'shape',@(res)tube_in_plate(ids{1},tt,ls{1},pt,37.5,res)), ...
struct('thickness',ts{2},'shape',@(res)tube_in_plate(ids{2},tt,ls{2},pt,37.5,res)), ...
struct('thickness',ts{3},'shape',@(res)tube_in_plate(ids{3},tt,ls{3},pt,37.5,res)), ...
struct('thickness',ts{4},'shape',@(res)tube_in_plate(ids{4},tt,ls{4},pt,37.5,res)), ...
struct('thickness',ts{5},'shape',@(res)tube_in_plate(ids{5},tt,ls{5},pt,37.5,res)), ...
struct('thickness',ts{6},'shape',@(res)tube_in_plate(ids{6},tt,ls{6},pt,37.5,res)), ...
{'space',1}, ...
struct('thickness',8,'shape',@(res)tube_in_plate(3,0.5,6,2,37.5,res)) ...
};

  % load pre‐computed quadrupole bender DA maps
  ld = load('MAT007_clean_DA/bender_segmented_data.mat', ...
            'transfer_maps','ref_states','z_coords');
  transfer_maps_bender = ld.transfer_maps;
  ref_states_bender    = ld.ref_states;
  z_coords_bender      = ld.z_coords;

  % build stack & run BEM once
  opts_stack = struct('z0',simParams.z0,'spacing',simParams.spacing,'res',simParams.res);
  [electrodes_full_stack, z_boundaries] = createElectrodeStack(electrodeData, opts_stack);

  bidx = find(cellfun(@(c)iscell(c)&&strcmpi(c{1},'space'), electrodeData),1);
  z_bend_start = z_boundaries(bidx) + (7.95 - (116 - 110.696));
  z_bend_end   = z_boundaries(bidx+1);
  z_coords_bender = z_coords_bender + z_bend_start;

  [qs, bemTable] = BEM_monopole(electrodes_full_stack, 1.0, false);
  bemParams.q_charge = daParams.q_charge;
  bemParams.m_mass   = daParams.m_mass;
  bemParams.E0       = daParams.E0;

  % build DA maps for bender
  sim_daOrder = daParams.daOrder;
  sim_nVars   = daParams.nVars;
  r_da  = DiffAlg.var(1,sim_daOrder,sim_nVars);
  pr_da = DiffAlg.var(2,sim_daOrder,sim_nVars);
  one_da= DiffAlg.one(sim_nVars,sim_daOrder);
  nseg = numel(transfer_maps_bender);
  S_maps_bender_segmented = cell(1,nseg);
  for i = 1:nseg
    M = transfer_maps_bender{i};
    X0 = ref_states_bender(:,i);
    X1 = ref_states_bender(:,i+1);
    Off = X1 - M*X0;
    S_maps_bender_segmented{i} = { ...
      r_da*M(1,1) + pr_da*M(1,2) + one_da*Off(1), ...
      r_da*M(2,1) + pr_da*M(2,2) + one_da*Off(2)  ...
    };
  end

  % ray‐trace parameters
  traceParams.num_r_rays = 16;
  traceParams.num_s_rays = 16;
  traceParams.max_r0_dev = 1.5;
  traceParams.max_s0_dev = 0.01;
  traceParams.center_r0  = 0.0;
  traceParams.center_s0  = 0.0;
  visParams.electrodes_full_stack = electrodes_full_stack;
  visParams.simParams = simParams;
  trace_eval_order = sim_daOrder;

  %% ----------------------
  % 2) SETUP SIMULATED ANNEALING + HYBRID fmincon
  %% ----------------------
  % bounds & start
  LB = -5*ones(8,1);
  UB =  0.5*ones(8,1);
  x0 = mean([LB,UB],2);

  % SA options
  saOpts = saoptimset( ...
    'Display','iter', 'MaxIter',1000, 'MaxFunEvals',50000, ...
    'TemperatureFcn',@temperaturefast, 'ReannealInterval',50, ...
    'PlotFcns',{@saplotbestx,@saplotf,@saplotstopping}, ...
    'HybridFcn',{@fmincon,optimoptions('fmincon',... 
       'Algorithm','sqp','Display','off','UseParallel',true, ...
       'TolFun',1e-8,'TolX',1e-8)} ...
  );

  fprintf('Starting Simulated Annealing...\n');
  [x_opt,fval,exitflag,output] = simulannealbnd(@objective, x0, LB, UB, saOpts);
  fprintf('SA complete: RMS = %.6g (exitflag=%d)\n', fval, exitflag);
  fprintf('Final voltages:\n el1=%.4f, se1=%.4f, se2=%.4f, dl1=%.4f, dl2=%.4f, dl3=%.4f, dl4=%.4f, ecv=%.4f\n', x_opt);

  if nargout<1, clear x_opt; end

  %% ----------------------
  % Nested OBJECTIVE: initial density-weighted
  %% ----------------------
  function rmsSize = objective(x)
    el1 = x(1); se1 = x(2);
    se2 = x(3); dl1 = x(4);
    dl2 = x(5); dl3 = x(6);
    dl4 = x(7); ecv = x(8);

    Vs = [1.0; 0.75; 0.50; 0.25; ...    % 1-4: uniform field tube
                0.0; el1; 0.0; ...            % 5-7: Einzel 1
                se1; 0.0; 0.0; 0.0; ...       % 8-11: steering electrode 1
                % No voltage for the bender 'space' (item 10)
                0.0; se2; 0.0;...             % 10-12: Einzel 2
                0.0; ...
                0.0; dl1; dl2; dl3; dl4; ecv; % Decelerator assembly
                -0.5; % endcap electrode
               ];
    bemParams.qVs = qs' * Vs;

    p1 = struct('z_start',dopriParams.z_start_overall,'z_end',z_bend_start);
    [z1,S1,~] = run_dopri87_da_integration(p1, daParams, bemParams, bemTable, integratorOpts);

    p2 = struct('z_start',z_bend_end,'z_end',dopriParams.z_end_overall);
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
                    false);
    % drawnow;

    % extract final radii & compute RMS spot size
    ri = squeeze(rpos(  1, :));
    rf = squeeze(rpos(end, :));
    gauss = @(r) exp(-(r / std(ri)).^2);
    % constants can be dropped.
    opval = abs(rf .* (ri .* gauss(ri)));
    rmsSize = std(opval);
    % fprintf("Objective function scoring: %.4g rms radius -> scored %.4g\n", ...
    %     std(rf), rmsSize)
  end
end
