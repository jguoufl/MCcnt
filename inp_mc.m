%%%%%%%% inputs for the simulation %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global q kBT hbar
global n t_step qsup Ntran Nstdy Degflag Balflag    % MC simulation parameters 
global Nkf Nxf                                      % E, k, and x grid

global Evec ga_ap ga_op ga_tot hw tau_IR     % the output of rate.m

%%%%%% constants
q=1.6e-19;
kBT=0.0259;
hbar=1.055e-34;

%%%%%%%%%%%%% inputs %%%%%%%%%%%%%%
n=25;                 % (n, 0) CNT
lamda_ap=10*1e-9;   % acoustic phonon mfp at high energy.
lamda_op=10000*1e-9;   % optical phonon mfp at high energy
hw=0.20;            % the optical phonon energy in eV
tau_IR=1*1e15;    % the spontaneous IR emission time
t_step=0.75*1e-15; % the time step in S
qsup=0.001;   % super-particle charge in e
Ntran=100; Nstdy=200;   % No. of transient and steady state time steps

%%%%% the grid for collecting carrier statistics
Degflag=1;  % 1 for treating Pauli exclusion, 0 for not treating
Balflag=0;  % ballitic transport flag
dE=10e-3;    % the energy bin in eV for collecting the current spectrum
Nxf=100;     % the number of position bins for collecting the distribution function, fd(x,k)
Nkf=100;    % number of bins for k
Np_max=0;   % the initial number of particles

