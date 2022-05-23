function simout = runVlmValidation( aircraft_name, gust, flaps )


%% create wing
switch aircraft_name
    case 'leisa'
        wing = wingCreate('wing_params_leisa_main',50,'spacing','like_chord','is_unsteady',true);
        fuselage = fuselageCreate('fuselage_params_leisa',4,20,'unsteady');
        
        % define current rigid body state
        alpha = deg2rad(1.49+0.5);
        V = 237.23;
        h = 10668;

        % cg position in wing coordinate system
        xyz_cg = [-19.5;0;0];
    case 'se2a'
        addPathTiGL('2.2.3')
        tixiHandle = tixiOpenDocumentTry( ... 
        which ( 'SE2A_AC_Design_MR_V2_BwdSweep_CPACS2_Turbulent.xml' ) );
        tiglHandle = tiglOpenCPACSConfigurationTry( tixiHandle );
        wing = wingCreateWithCPACS( tiglHandle, 1, 50, 'spacing', 'like_chord', 'airfoil_method', 'analytic', 'is_unsteady', true );
        axis_reversed = [ -1; 1; -1 ];
        fuselage = fuselageCreateWithCpacs( tiglHandle, 'Fuse', axis_reversed, 20, 'unsteady' );
        
        alpha = deg2rad(-0.5);
        V = 240;
        h = 6000;
        xyz_cg = [-19.5;0;0];
        
end

omega = [0;0;0];
beta = 0;
V_Kb_dt = zeros(3,1);
omega_dt = zeros(3,1);


lift2bm = wingGetBendingMoment( wing.geometry, 0.05 );


% define actuator states
actuators_main_pos = zeros( size( flaps.magn ) );
actuators_main_rate = zeros( size( flaps.magn ) );


%% structure state
% structure_state = zeros(2*size(structure_red.K,1),1);
structure_state = zeros(2*1,1);

structure_accel = zeros(1,1);

%% aerodynamic state
unst_aero_x = zeros(size(wing.state.aero.unsteady.x));
unst_aero_X = zeros(size(wing.state.aero.unsteady.X));
unst_aero_z = zeros(size(wing.state.aero.unsteady.z));
unst_aero_z2 = zeros(size(wing.state.aero.unsteady.z2));
alpha_ind = zeros(size(wing.state.aero.circulation.alpha_ind));
alpha_unst = zeros(size(fuselage.state.aero.unsteady.alpha));
beta_unst = zeros(size(fuselage.state.aero.unsteady.beta));
gamma_filt = zeros( 5, size(wing.state.aero.circulation.alpha_ind,2) );
wing.state.aero.circulation.gamma_filt = gamma_filt;

%% atmosphere
atmosphereStruct = isAtmosphere(h);

%% sample time
dt = 0.001;

%% init state
t = 0;

%% gust

f = 2*pi/(gust.lambda/V);
U = zeros( 3, wing.n_panel );
U_dt = zeros( 3, wing.n_panel );
U3 = zeros( 1, wing.n_panel );
U3_dt = zeros( 1, wing.n_panel );
V_Wb = zeros( 3, fuselage.n_segments + 1 );
V_Wb_dt = zeros( 3, fuselage.n_segments + 1 );
V3_Wb = zeros( 1, fuselage.n_segments + 1 );
V3_Wb_dt = zeros( 1, fuselage.n_segments + 1 );
actuators_main_pos1 = zeros( size( wing.state.actuators.pos ) );
actuators_main_rate1 = zeros( size( wing.state.actuators.pos ) );

%% init plot values
num_samples = 1500;
simout.c_L = repmat( wing.state.aero.coeff_loc.c_XYZ_b(3,:), num_samples, 1 );
simout.Delta_alpha = repmat( wing.state.aero.coeff_loc.c_XYZ_b(3,:), num_samples, 1 );
simout.C_XYZ_b = repmat( wing.state.aero.coeff_glob.C_XYZ_b, 1, num_samples );
simout.C_bm = zeros( num_samples, size(lift2bm,1) );
simout.C_XYZ_fuse_b = repmat( wing.state.aero.coeff_glob.C_XYZ_b, 1, num_samples );
simout.alpha_unst = zeros( length(fuselage.state.aero.unsteady.alpha), num_samples );
simout.num_iter = zeros(1,num_samples);


%% compute wing state

delta_t_disp = 0.1;
t_disp = 0.1;

disp('VLM simulation startet.')
tic;
for k = 1:num_samples
    
    if t >= t_disp
        disp(['t=',num2str(t),'s']);
        t_disp = t_disp + delta_t_disp;
    end
    
    t0_shift = gust.t0 - ( wing.geometry.ctrl_pt.pos(1,:) + wing.geometry.origin(1)*0 )/V;
    is_no_gust = t < t0_shift | t > t0_shift + 2*pi/f;    
    U3( is_no_gust ) = 0;
    U3_dt( is_no_gust ) = 0;
    U3( ~is_no_gust ) = -gust.U_ds/2* (1-cos(f*(t-t0_shift(~is_no_gust))));
    U3_dt( ~is_no_gust ) = -gust.U_ds/2*f*sin(f*(t-t0_shift(~is_no_gust)));
    U(3,:) = U3;
    U_dt(3,:) = U3_dt;

    t0_shift_flap = gust.t0;
    if t < t0_shift_flap || t> t0_shift_flap + 2*pi/f
        actuators_main_pos(:) = 0;
        actuators_main_rate(:) = 0;
    else
        actuators_main_pos(:) = -flaps.magn/2 .* (1-cos(flaps.freq*(t-t0_shift_flap)));
        actuators_main_rate(:) = -flaps.magn/2 .* flaps.freq .* sin(flaps.freq.*(t-t0_shift_flap))*0;
    end
    
%     actuators_main_pos(:) = flaps.magn .* sin( 2*pi * flaps.freq * t );
%     actuators_main_rate(:) = 2*pi * flaps.freq .* flaps.magn .* cos( 2*pi * flaps.freq * t );

    wing_out = wingSetState(wing, alpha, beta, ...
        V, omega, actuators_main_pos, actuators_main_rate, xyz_cg, ...
        'V_Kb_dt', V_Kb_dt, 'omega_dt', omega_dt, ...
        'atmosphere', atmosphereStruct, 'wind', U, U_dt, ...
        'structure_pos', structure_state(1:end/2), 'structure_vel', structure_state(end/2+1:end), ...
        'structure_accel', structure_accel, ...
        'unst_airfoil_state', unst_aero_x, 'dyn_stall_state', unst_aero_X, ...
        'unst_flap_state', unst_aero_z, 'unst_act2_state', unst_aero_z2, ...
        'alpha_ind', alpha_ind );
    
%     wing = wingSetState( wing, alpha, beta, V, omega, actuators_pos, actuators_rate, [0;0;0], ...
%         'atmosphere', atmosphereStruct );
    
    % bending moment
    [ XYZ_i_b ]     = wingGetLocalForce( wing_out );
    lift = -XYZ_i_b(3,:);
    bm = lift2bm * lift';
    C_bm = bm/(0.5*wing_out.state.external.atmosphere.rho*V^2*wing_out.params.S*wing_out.params.b/2);

    % time integration (Euler forward)
    unst_aero_x = wing_out.state.aero.unsteady.x + wing_out.state.aero.unsteady.x_dt * dt;
    unst_aero_X = wing_out.state.aero.unsteady.X + wing_out.state.aero.unsteady.X_dt * dt;
    unst_aero_z = wing_out.state.aero.unsteady.z + wing_out.state.aero.unsteady.z_dt * dt;
    unst_aero_z2 = wing_out.state.aero.unsteady.z2 + wing_out.state.aero.unsteady.z2_dt * dt;
    Ma = min(0.95,V / wing_out.state.external.atmosphere.a);
    beta_Ma = sqrt( 1 - (Ma*cosd(30))^2 );
    c = wing_out.params.S/wing_out.params.b;
    b12 = 0.45;
    
%     T_downwash = ( 1./wing.state.geometry.ctrl_pt.c + 0*(wing.state.geometry.ctrl_pt.pos(1,:)+wing.state.geometry.ctrl_pt.c/2) - min( (wing.state.geometry.ctrl_pt.pos(1,:)+wing.state.geometry.ctrl_pt.c/2) ) ) / V;
    
%     T_downwash = 0.01*1 ./ ( (2*V./c^3)*beta_Ma^2*b12 ) .* wing.state.geometry.ctrl_pt.c.^2;
%     T_downwash = wing.state.geometry.ctrl_pt.c.^2 / V;
    T_downwash = 2*c/V;
%     alpha_ind = alpha_ind + 1./T_downwash .* ( wing.state.aero.circulation.alpha_ind - alpha_ind ) * dt;
    alpha_ind = wing_out.state.aero.circulation.alpha_ind;
    gamma_filt = gamma_filt + 1./T_downwash .* ( wing_out.state.aero.circulation.gamma - gamma_filt ) * dt;
    wing_out.state.aero.circulation.gamma_filt = gamma_filt;
    
    t = t + dt;
    
    % output
    c_XYZ_a = dcmBaFromAeroAngles(alpha,beta)' * wing_out.state.aero.coeff_loc.c_XYZ_b;
    simout.c_L(k,:) = -c_XYZ_a(3,:);
    simout.C_XYZ_b(:,k) = wing_out.state.aero.coeff_glob.C_XYZ_b;
    simout.num_iter(k) = wing_out.state.aero.circulation.num_iter;
    simout.Delta_alpha(k,:) = wing_out.state.aero.circulation.Delta_alpha;
    
    simout.C_bm(k,:) = C_bm;
    
    % fuselage
    t0_shift_fuse = gust.t0 - ( fuselage.geometry.border_pos(1,:) - wing_out.geometry.origin(1) )/V;
    is_no_gust = t < t0_shift_fuse | t > t0_shift_fuse + 2*pi/f;    
    V3_Wb( is_no_gust ) = 0;
    V3_Wb_dt( is_no_gust ) = 0;
    V3_Wb( ~is_no_gust ) = gust.U_ds/2* (1-cos(f*(t-t0_shift(~is_no_gust))));
    V3_Wb_dt( ~is_no_gust ) = gust.U_ds/2*f*sin(f*(t-t0_shift(~is_no_gust)));
    V_Wb(3,:) = V3_Wb;
    V_Wb_dt(3,:) = V3_Wb_dt;
    
    fuselage = fuselageSetState( fuselage,alpha,beta,V,omega,xyz_cg,...
        'wind', V_Wb, V_Wb_dt, 'unsteady', alpha_unst, beta_unst );
%     fuselage = fuselageSetState( fuselage,alpha,beta,V,omega,xyz_cg,...
%         'wind', V_Wb, V_Wb_dt );
    
    alpha_unst = alpha_unst + fuselage.state.aero.unsteady.alpha_dt * dt;
    
    simout.alpha_unst(:,k) = fuselage.state.aero.unsteady.alpha';

    simout.C_XYZ_fuse_b(:,k) = fuselage.state.aero.R_Ab / (0.5*wing_out.state.external.atmosphere.rho*V^2*wing_out.params.S);
    
end
elapsed_time = toc;
disp(['VLM simulation finished after ',num2str(elapsed_time),' seconds.'])

simout.time = 0:dt:t;
simout.wing = wing_out;

end