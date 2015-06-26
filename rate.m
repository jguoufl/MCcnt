%%%%%%% scattering rates vs. energies %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Evec, ga_ap, ga_op, ga_ii, ga_tot]=rate(lamda_ap,lamda_op,tau_IR,Sii,hw,Egh1)
%% input: tau_IR: constant IR recombination lifetime
%%        Sii: in /s, the rate parameter for impact ionization
%%        n: integer, for (n,0) CNT index
%% output: Evec: the kinetic energy vector
%%         ga*: scattering rates

Eth=2*Egh1;     % the threshold for impact ionization
Ekmax=1;     % eV,the maximum kinetic energy considered in scattering and carrier statistics
vF=9.8e5;               % m/s, Fermi velocity
Ek=inline('Egh*(sqrt(1+(k./k0).^2)-1)','k','Egh','k0');
vE=inline('vF*sqrt(1-(1./(1+E./Egh)).^2)','E','Egh','vF');
kE=inline('k0*sqrt(((E./Egh)+1).^2-1)','E','Egh','k0');   
kv=inline('k0*v./sqrt(vF^2-v.^2)','v','vF','k0');   

brd=2.5e-3;  % the broadening of DOS due to phonon perturbation in eV.
NE=1001;   % number of energy grids for computing DOS
Evec=linspace(0,Ekmax,NE);  % the energy grid for the scattering rate
Edos=linspace(-Ekmax,Ekmax,2*NE-1);   % the energy grid for DOS
dE_ga=Evec(2)-Evec(1);
lim_up=Ekmax+Egh1+5*brd;   % upper limit of the integral
ff=inline('exp(-(E+Egh-sqrt(x^2+Egh^2)).^2./(2*brd^2))','x','E','Egh','brd');
Dn=(1/(sqrt(2*pi)*brd))*quadv(ff,0,lim_up,1e-5,[],Edos,Egh1,brd);

%ga_ap=(1/lamda_ap)*vE(Evec,Egh1,vF);  % hypothetical elatic phonon with constant ga;
ga_ap=(vF/lamda_ap)*interp1(Edos,Dn,Evec);  % acoustic phonon scattering rate, 300nm mfp
ga_op=(vF/lamda_op)*interp1(Edos,Dn,Evec-hw);   % zb1 phonon scattering rate, 15nm mfp
ga_ii=Sii*((Evec-Eth)./Eth).^2.*heaviside(Evec-Eth);     % Impact ionization rate
ga_tot=ga_ap+ga_op+ga_ii+1/tau_IR;  % total scattering rate
ga_max=max(ga_tot);
%%%%%%%%%% end of computing scattering rate vs. energy %%%%%%%