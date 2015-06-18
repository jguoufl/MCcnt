%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non self-consistent Monte-Carlo simulation of CNTFETs

clear all
close all

inp_mc;
[Evec, ga_ap, ga_op, ga_tot]=rate(lamda_ap,lamda_op, tau_IR,hw,n);
load feton  % load the first subband profile for non-self-consistent simulation
Vd_bias=0.6;
[Ne_bias,Id, jS, jD, Emesh]=mc(XI,E1,Vd_bias);
draw