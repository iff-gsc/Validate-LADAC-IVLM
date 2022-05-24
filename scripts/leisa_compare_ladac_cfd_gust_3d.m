% comparison of VLM with RANS simulation during 1-cos gust encounter for
% LEISA aircraft

%% add folders to path
addPathIvlmValidation();


%% load RANS simulation data
cfd_3d = cfdJunaidImportWingAll('GustOnly/wing');

%% run VLM simulation

% peak 1-cos gust velocity, in m/s
gust.U_ds = 14.49;
% gust length, in m
gust.lambda = 50;
% time before gust hits aircraft nose, in s
gust.t0 = 0.0421;

% aircraft velocity, in m/s
state.V = 237.23;
% angle of attack, in rad
state.alpha = deg2rad(1.49+0.5);
% altitude, in m
state.h = 10668;

% flap settings
flaps.freq = zeros(1,5);
flaps.magn = zeros(1,5);

vlmout = runVlmValidation( 'leisa', gust, state, flaps );

%% use the same sampling for VLM and RANS

c_L_VLM_interp = zeros( length(cfd_3d.time), vlmout.wing.n_panel );
for i = 1:vlmout.wing.n_panel
    c_L_VLM_interp(:,i) = interp1( vlmout.time, vlmout.c_L(:,i), cfd_3d.time )';
end

eta_interp_idx = vlmout.wing.geometry.ctrl_pt.pos(2,:)/vlmout.wing.params.b*2 >= min( cfd_3d.eta(1,:) );
eta_interp = vlmout.wing.geometry.ctrl_pt.pos(2,eta_interp_idx)/vlmout.wing.params.b*2;
c_L_RANS_interp = zeros( length(cfd_3d.time), length( eta_interp ) );
for i = 1:size(cfd_3d.eta,1)
    c_L_RANS_interp(i,:) = interp1( cfd_3d.eta(i,:), cfd_3d.cl(i,:), eta_interp );
end

[X1,Y1] = meshgrid(cfd_3d.time,eta_interp);

error_mat = c_L_RANS_interp - c_L_VLM_interp( :, eta_interp_idx );

%% compare steady local c_L

figure
plot(cfd_3d.Y(20,:)/cfd_3d.Y(20,end),cfd_3d.cl(2,:))
hold on
plot(vlmout.wing.geometry.ctrl_pt.pos(2,:)/vlmout.wing.params.b*2, -vlmout.wing.state.aero.coeff_loc.c_XYZ_b(3,:) )
grid on
xlim([0,1])
ylim([0,0.9])
xlabel('Dimensionless span coordinate')
ylabel('Local lift coefficient')
legend('RANS (TAU Code)', 'VLM (LADAC)')
title('Steady local lift coefficient comparison')

%% compare local c_L

time_idx_end = 200;

figure
[X,Y] = meshgrid(cfd_3d.time(1:time_idx_end),eta_interp);
surf(X,Y,c_L_RANS_interp(1:time_idx_end,:)','FaceAlpha',0.3)
hold on
surf(X,Y,c_L_VLM_interp(1:time_idx_end,eta_interp_idx)','FaceAlpha',0.3,'EdgeColor','r')
grid on
ylabel('Dimensionless span, -')
xlabel('Time, s')
zlabel('Local lift coefficient, -')
legend('RANS (TAU Code)', 'VLM (LADAC)')
title('Local lift coefficient comparison during 1-cos gust')

%% show absolute local c_L error

figure
surf(X1,Y1,error_mat','FaceAlpha',0.3)
xlabel('Time, s')
ylabel('Dimensionless span')
zlabel('Absolute local lift coefficient error')
title('Local lift coefficient error during 1-cos gust')

%% compare WRBM coefficient

figure
c_interp = vlmout.wing.state.geometry.ctrl_pt.c(eta_interp_idx);
% c_interp(:) = 1;
b_interp = diff( [ eta_interp(1), eta_interp(1:end-1)+diff(eta_interp)/2, eta_interp(end) ] )*vlmout.wing.params.b/2;
plot( cfd_3d.time, c_L_RANS_interp * ( eta_interp .* c_interp .* b_interp )' / vlmout.wing.params.S )
hold on
plot( cfd_3d.time, c_L_VLM_interp(:,eta_interp_idx) * ( eta_interp .* c_interp .* b_interp )' / vlmout.wing.params.S )
grid on
xlabel('Time, s')
ylabel('Wing root bending moment coefficient')
legend('RANS (TAU Code)','VLM (LADAC)')
title('Comparison of wing root bending moment coefficient during 1-cos gust')
