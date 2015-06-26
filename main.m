%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non self-consistent Monte-Carlo simulation of CNTFETs

clear all
close all

%% input of  the program
inp_mc;
%% compute the scattering rate
[Evec, ga_ap, ga_op,ga_ii, ga_tot]=rate(lamda_ap,lamda_op,tau_IR,Sii0,hw,Egh1);
%%% load E1(XI) with a bias of Bd_bias
load profile2  % load the first subband profile for non-self-consistent simulation
%%% MC simulation
[Ektp xp vp Ne_bias,Id, jS, jD, Emesh]=mc(XI,E1,Vd_bias);
%% visualization
draw