%% init

% add folders to path
addPathIvlmValidation();


%% load RANS simulation data

% cfd_3d = cfdKhalidImportWingAll('se2a-data/distributions');

%% run VLM simulation

gust.U_ds = 14.49;
gust.lambda = 50;
gust.t0 = 0.8;

flaps.freq = zeros(1,5);
flaps.magn = zeros(1,5);

vlmout = runVlmValidation( 'se2a', gust, flaps );

%% use the same sampling for VLM and RANS

idx = (vlmout.time + 0.1) > gust.t0;

c_L_VLM_interp = zeros( length(cfd_3d.time), vlmout.wing.n_panel );
for i = 1:vlmout.wing.n_panel
    c_L_VLM_interp(:,i) = interp1( vlmout.time(idx)-gust.t0+0.1, vlmout.c_L(idx,i), cfd_3d.time )';
end

eta_interp_idx = vlmout.wing.geometry.ctrl_pt.pos(2,:)/vlmout.wing.params.b*2 >= min( cfd_3d.eta(1,:) );
eta_interp = vlmout.wing.geometry.ctrl_pt.pos(2,eta_interp_idx)/vlmout.wing.params.b*2;
c_L_RANS_interp = zeros( length(cfd_3d.time), length( eta_interp ) );

V = 240.463;
h = 6000;
rho = 0.6597;
S = 79.34688;

Y_diff = [diff(cfd_3d.Y(1,1:2)),diff(cfd_3d.Y(1,:))];
S_local = Y_diff .* interp1( vlmout.wing.geometry.ctrl_pt.pos(2,end/2:end), vlmout.wing.geometry.ctrl_pt.c(end/2:end), cfd_3d.Y(1,:), 'linear', 'extrap' );
for i = 1:size(cfd_3d.eta,1)
    c_L_RANS_interp(i,:) = interp1( cfd_3d.eta(i,:), cfd_3d.fz(i,:).*Y_diff./(0.5*rho*V^2.*S_local), eta_interp );
end
% c_L_RANS_interp = c_L_RANS_interp ./ ( vlmout.wing.geometry.ctrl_pt.c(end/2+1:end) .* diff(vlmout.wing.geometry.vortex.pos(2,end/2:end)) );

[X1,Y1] = meshgrid(cfd_3d.time,eta_interp);

error_mat = c_L_RANS_interp - c_L_VLM_interp( :, eta_interp_idx );

%% compare steady local c_L

figure
plot(eta_interp,c_L_RANS_interp(1,:))
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

time_idx_end = length(cfd_3d.time);

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
