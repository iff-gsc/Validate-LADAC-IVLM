function simout = runVlmValidation( aircraft_name, gust, state, flaps )
% runVlmValidation simulate VLM gust encounters and flaps deflections for a
% fuselage wing combination that is fixed (no rigid body motion and no
% structural deflections) and compare the results with RANS simulation data

disp('VLM simulation startet.')
tic;

%% simulation parameters
% sample time
dt = 0.001;
% display simulation progress in command window each delta_t_disp seconds
delta_t_disp = 0.1;
% simulation stop time, in s
t_end = 1.5;

%% init wing and fuselage
switch aircraft_name
    case 'leisa'
        wing = wingCreate('wing_params_leisa_main',50,'spacing','like_chord','is_unsteady',true);
        wing_static = wingCreate('wing_params_leisa_main',50,'spacing','like_chord','is_unsteady',false);
        fuselage = fuselageCreate('fuselage_params_leisa',4,20,'unsteady');
       
        % cg position in wing coordinate system
        xyz_cg = [-19.5;0;0];
    case 'se2a'
        addPathTiGL('2.2.3')
        tixiHandle = tixiOpenDocumentTry( ... 
        which ( 'SE2A_AC_Design_MR_V2_BwdSweep_CPACS2_Turbulent.xml' ) );
        tiglHandle = tiglOpenCPACSConfigurationTry( tixiHandle );
        wing = wingCreateWithCPACS( tiglHandle, 1, 50, 'spacing', 'like_chord', 'airfoil_method', 'analytic', 'is_unsteady', true );
        wing_static = wingCreateWithCPACS( tiglHandle, 1, 50, 'spacing', 'like_chord', 'airfoil_method', 'analytic', 'is_unsteady', false );
        axis_reversed = [ -1; 1; -1 ];
        fuselage = fuselageCreateWithCpacs( tiglHandle, 'Fuse', axis_reversed, 20, 'unsteady' );
        
        xyz_cg = [-19.5;0;0];
        
end

omega = [0;0;0];
beta = 0;

lift2bm = wingGetBendingMoment( wing.geometry, 0.05 );

% define actuator states
actuators_main_pos = zeros( size( flaps.magn ) );
actuators_main_rate = zeros( size( flaps.magn ) );


%% atmosphere
atmosphere_struct = isAtmosphere(state.h);

%% gust
f = 2*pi/(gust.lambda/state.V);
U = zeros( 3, wing.n_panel );
U_dt = zeros( 3, wing.n_panel );
U3 = zeros( 1, wing.n_panel );
U3_dt = zeros( 1, wing.n_panel );
V_Wb = zeros( 3, fuselage.n_segments + 1 );
V_Wb_dt = zeros( 3, fuselage.n_segments + 1 );
V3_Wb = zeros( 1, fuselage.n_segments + 1 );
V3_Wb_dt = zeros( 1, fuselage.n_segments + 1 );

%% init output values
num_samples = floor(t_end/dt);
simout.c_L = repmat( wing.state.aero.coeff_loc.c_XYZ_b(3,:), num_samples, 1 );
simout.Delta_alpha = repmat( wing.state.aero.coeff_loc.c_XYZ_b(3,:), num_samples, 1 );
simout.C_XYZ_b = repmat( wing.state.aero.coeff_glob.C_XYZ_b, 1, num_samples );
simout.C_bm = zeros( num_samples, size(lift2bm,1) );
simout.C_XYZ_fuse_b = repmat( wing.state.aero.coeff_glob.C_XYZ_b, 1, num_samples );
simout.alpha_unst = zeros( length(fuselage.state.aero.unsteady.alpha), num_samples );
simout.num_iter = zeros(1,num_samples);

%% compute initial steady-state wing states
wing_static = wingSetState(wing_static, state.alpha, 0, state.V, ...
    zeros(3,1), 0, 0, zeros(3,1), 'atmosphere', atmosphere_struct );

%% init aerodynamic state
% steady state wing
[ unst_aero_x, unst_aero_X, unst_aero_z, unst_aero_z2] = ...
    wingStateGetUnstAeroTrimOp( wing_static.state, ...
    wing_static.airfoil, wing_static.config );
alpha_ind = wing_static.state.aero.circulation.alpha_ind;

% fuselage
alpha_unst = zeros(size(fuselage.state.aero.unsteady.alpha));
beta_unst = zeros(size(fuselage.state.aero.unsteady.beta));

%% unsteady simulation
% init time
t = 0;
% init display time
t_disp = 0.1;
% simulate
for k = 1:num_samples
    
    % display simulation progress
    if t >= t_disp
        disp(['t=',num2str(t),'s']);
        t_disp = t_disp + delta_t_disp;
    end
    
    % wing: set wind velocity in case there is a gust
    t0_shift = gust.t0 - ( wing.geometry.ctrl_pt.pos(1,:) + wing.geometry.origin(1) )/state.V;
    is_no_gust = t < t0_shift | t > t0_shift + 2*pi/f;    
    U3( is_no_gust ) = 0;
    U3_dt( is_no_gust ) = 0;
    U3( ~is_no_gust ) = -gust.U_ds/2* (1-cos(f*(t-t0_shift(~is_no_gust))));
    U3_dt( ~is_no_gust ) = -gust.U_ds/2*f*sin(f*(t-t0_shift(~is_no_gust)));
    U(3,:) = U3;
    U_dt(3,:) = U3_dt;
    
    % fuselage: set wind velocity in case there is a gust
    t0_shift_fuse = gust.t0 - fuselage.geometry.border_pos(1,:)/state.V;
    is_no_gust = t < t0_shift_fuse | t > t0_shift_fuse + 2*pi/f;    
    V3_Wb( is_no_gust ) = 0;
    V3_Wb_dt( is_no_gust ) = 0;
    V3_Wb( ~is_no_gust ) = gust.U_ds/2* (1-cos(f*(t-t0_shift(~is_no_gust))));
    V3_Wb_dt( ~is_no_gust ) = gust.U_ds/2*f*sin(f*(t-t0_shift(~is_no_gust)));
    V_Wb(3,:) = V3_Wb;
    V_Wb_dt(3,:) = V3_Wb_dt;

    % wing: set flap deflection if it was defined
    t0_shift_flap = gust.t0;
    if t < t0_shift_flap || t> t0_shift_flap + 2*pi/f
        actuators_main_pos(:) = 0;
        actuators_main_rate(:) = 0;
    else
        actuators_main_pos(:) = -flaps.magn/2 .* (1-cos(flaps.freq*(t-t0_shift_flap)));
        actuators_main_rate(:) = -flaps.magn/2 .* flaps.freq .* sin(flaps.freq.*(t-t0_shift_flap))*0;
    end

    % wing: compute unsteady aerodynamics
    wing_out = wingSetState(wing, state.alpha, beta, ...
        state.V, omega, actuators_main_pos, actuators_main_rate, xyz_cg, ...
        'atmosphere', atmosphere_struct, 'wind', U, U_dt, ...
        'unst_airfoil_state', unst_aero_x, 'dyn_stall_state', unst_aero_X, ...
        'unst_flap_state', unst_aero_z, 'unst_act2_state', unst_aero_z2, ...
        'alpha_ind', alpha_ind );
    
    % fuselage: compute unsteady aerodynamics
    fuselage = fuselageSetState( fuselage,state.alpha,beta,state.V,omega,xyz_cg,...
        'wind', V_Wb, V_Wb_dt, 'unsteady', alpha_unst, beta_unst );
    
    % wing: get bending moment coefficient
    [ XYZ_i_b ]     = wingGetLocalForce( wing_out );
    lift = -XYZ_i_b(3,:);
    bm = lift2bm * lift';
    C_bm = bm/(0.5*wing_out.state.external.atmosphere.rho*state.V^2*wing_out.params.S*wing_out.params.b/2);

    % wing: 1st order delay downwash
    c = wing_out.params.S/wing_out.params.b;
    T_downwash = 2*c/state.V;
    alpha_ind_dt = 1./T_downwash .* ( wing_out.state.aero.circulation.alpha_ind - alpha_ind );
    
    % wing: time integration (Euler forward)
    unst_aero_x = wing_out.state.aero.unsteady.x + wing_out.state.aero.unsteady.x_dt * dt;
    unst_aero_X = wing_out.state.aero.unsteady.X + wing_out.state.aero.unsteady.X_dt * dt;
    unst_aero_z = wing_out.state.aero.unsteady.z + wing_out.state.aero.unsteady.z_dt * dt;
    unst_aero_z2 = wing_out.state.aero.unsteady.z2 + wing_out.state.aero.unsteady.z2_dt * dt;
    alpha_ind = alpha_ind + alpha_ind_dt * dt;
    
    % fuselage: time integration (Euler forward)
    alpha_unst = alpha_unst + fuselage.state.aero.unsteady.alpha_dt * dt;
    
    % increment time
    t = t + dt;
    
    % set outputs
    c_XYZ_a = dcmBaFromAeroAngles(state.alpha,beta)' * wing_out.state.aero.coeff_loc.c_XYZ_b;
    simout.c_L(k,:) = -c_XYZ_a(3,:);
    simout.C_XYZ_b(:,k) = wing_out.state.aero.coeff_glob.C_XYZ_b;
    simout.num_iter(k) = wing_out.state.aero.circulation.num_iter;
    simout.Delta_alpha(k,:) = wing_out.state.aero.circulation.Delta_alpha;
    simout.C_bm(k,:) = C_bm;
    simout.alpha_unst(:,k) = fuselage.state.aero.unsteady.alpha';
    simout.C_XYZ_fuse_b(:,k) = fuselage.state.aero.R_Ab / ...
        (0.5*wing_out.state.external.atmosphere.rho*state.V^2*wing_out.params.S);
    
end

simout.time = 0:dt:t;
simout.wing = wing_out;

elapsed_time = toc;
disp(['VLM simulation finished after ',num2str(elapsed_time),' seconds.'])

end