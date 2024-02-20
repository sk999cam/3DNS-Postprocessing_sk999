%number of grid points in cube
ni = 584;
nj = 415;
nk = 96;
%shortest box edge length
minL =  0.009465;
%# of Fourier modes
M = 10000;
%Turbulence intensity
Tu = 0.05;
%Inlet/reference velocity
vIn = 42.9;
%vTu = sum(vrel.*rad_in)/sum(rad_in)
%Lengthscale of largest eddies 
L = -1;
%Dissipation
eps = 1.5e4;

%Calculate kniematic viscosity:
rgas = 2.871428571428571e+02;
Tor = 299.1747;
por = 1.1590e+05;
muref = 2.7360e-05;
T = Tor-0.4/(2*1.4*rgas)*vIn^2;
ro = por/(rgas*T)*(T/Tor)^(1.4/0.4);
nu = 1/ro*muref*((273+110.4)/(T+110.4))*(T/273)^1.5;


