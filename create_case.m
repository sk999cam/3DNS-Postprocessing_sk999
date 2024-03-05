% Soham Kar, University of Cambridge, 14 Feb 2024

% This script takes input from the used and creates a flat plate case to be
% run on 3DNS

casename = "cases/04032024_tripped_slow";

%%%  Meshing
L_ref = 1;

%create mesh
blk = mesh_channel_free(0.15/L_ref, 0.15/L_ref, 256, 384, 1280, 10);

%plot and check mesh
%mesh_analysis(blk);

% set solver inputs


%% boundary conditions
bcs.Toin = 300;
bcs.Poin = 200000;
bcs.pexit = 198600;
bcs.vin = 34.7;
bcs.alpha = 0.0;
bcs.gamma = 0.0;
bcs.aturb = 0.0;
bcs.lturb = 0.0;
bcs.readprof = 0;
bcs.twall = -1;
bcs.cax = 1;
bcs.theta = 2e-3/L_ref;

% freestream bc points
bcs.nfsp = 5;
bcs.k = 0.0e-6;
bcs.Lref = L_ref;

%dont understand these
bcs.nturb = 500;
bcs.iradprof = 0; %radial profile not relavent
bcs.g_z = 0; %spanwise body force set to 0

%% gas properties
gas.gamma = 1.4;
gas.cp = 1005.0;
gas.mu_ref = 1.0e-3;
gas.mu_tref = 273.0;
gas.mu_cref = 110.4;
gas.pr = 0.72;

%% solver inputs
solver.niter = 5000;
solver.nwrite = 1000;
solver.ncut = 1000000;
solver.cfl = 0.2;
solver.sigma = 0.05;
solver.irestart = 0;
solver.span = 0.2;
solver.fexpan = 1.0;
solver.nk = 1;
solver.npp = 1;
solver.istats = 2;
solver.version = ("gpu");
solver.ilam = 1;
solver.iblrec = 0;
solver.freebuff = 10;

% dont understand these inputs
solver.ifsplit = 1; %skew split differencing set to 1
solver.ifsat = 0; %
solver.ifLES = 0; % set to 0
solver.istability = 0; % set to 0

%% trip properties
trip.ifLEtrip = 0;
trip.x1 = 0;
trip.y1 = 0;
trip.x2 = 0;
trip.y2 = 0;
trip.tripscale = 0;
trip.kspace = 0;



%% write case files
write_input_files(casename, blk, bcs, gas, solver, trip)
write_case(casename, blk, gas, solver);

%channel = DNS_channel(casename);
%flo = channel.init_shock_flow(0.8, 0.5, 0.2, 2e-3);
%channel.instFlow = flo;
%channel.instFlow.writeFlow




