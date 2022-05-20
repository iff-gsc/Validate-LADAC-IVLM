% compute lift distribution of LEISA with Tornado

%% add folders to path
addPathIvlmValidation();

%% define wing
nx = 3;
% to do: wing camber is not considered!
geo = tornadoCreateGeo( 'wing_params_leisa_main',zeros(3,1),zeros(3,1),100,nx,3,zeros(1,5),true(1,5));

%% define wing state
state=tornadoCreateState(deg2rad(1.5),0,237.3*0.8,zeros(3,1),10667,1);
lattice_type = 1;
[lattice,ref]=fLattice_setup2(geo,state,lattice_type);

%% plot wing geometry
geometryplot(lattice,geo,ref);

%% run VLM
clear results
results.a=0;
results=solver9(results,state,geo,lattice,ref);
[results]=coeff_create3(results,lattice,state,ref,geo);

%% plot spanwise lift distribution
figure
plot( results.ystation/max(results.ystation), results.CL_local )
grid on
xlabel('Dimensionless span, -')
ylabel('Local lift coefficient, -')
ylim([0,inf])

