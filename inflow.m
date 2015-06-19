function [Eps_in Epd_in Ne_s Ne_d xp vp Ektp]=inflow(XI,E1,Emesh,mu_s,mu_d,Ne_s,Ne_d,xp,vp,Ektp)
%% Input: Emesh: energy mesh in Monte carlo
%%        E1: the first conduction subband edge
%%        Ne_s: source injection, [] for the 1st step, residue from the previous step otherwise
%%        Ne_d: drain injection
global kBT q hbar
global t_step qsup Egh1 vF

Lch=max(XI); % the total channel length
%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ek=inline('Egh*(sqrt(1+(k./k0).^2)-1)','k','Egh','k0');
vE=inline('vF*sqrt(1-(1./(1+E./Egh)).^2)','E','Egh','vF');
kE=inline('k0*sqrt(((E./Egh)+1).^2-1)','E','Egh','k0');   
kv=inline('k0*v./sqrt(vF^2-v.^2)','v','vF','k0'); 

Etops=E1(1);
Etopd=E1(end);
dE=Emesh(2)-Emesh(1);
Emesh_injs=Emesh(find(Emesh>Etops));                % inject above E1(1)
Emesh_injd=Emesh(find(Emesh>Etopd));                % inject above E1(1)
if length(Ne_s)==0  % intialization for the 1st time step
    Ne_s=zeros(1,length(Emesh_injs));                    % initialize the injection carrier spectrum
    Ne_d=zeros(1,length(Emesh_injd));
end


Ninjs=(t_step*4*kBT*q/(2*pi*hbar)*log(1+exp((mu_s-Etops)/kBT))/qsup);         % no. of S-injected superpaticles in each step
Ninjd=(t_step*4*kBT*q/(2*pi*hbar)*log(1+exp((mu_d-Etopd)/kBT))/qsup);   % no. of D-injected superpaticles in each step
Ninjt=Ninjs+Ninjd;                                                              % total no. of superpaticles in each step

infs_x=1e-10;                                       % in m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% inject carriers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %clear xp_inj Ektp_inj vp_inj 
    Eps_in=[]; Epd_in=[]; % initialization
    %%%%%%%%%%%%%%% the source injection, deterministic injection
    Ne_s=Ne_s+Ninjs/sum(1./(1+exp((Emesh_injs-mu_s)./kBT)))*(1./(1+exp((Emesh_injs-mu_s)./kBT)));   % the source injection spectrum
    mInj=find(Ne_s>=1);                             % the energy indices with carrier injcetion for Emesh_injs
    Ntrs=sum(floor(Ne_s(mInj)));                    % number of electrons injected from the source
    ind_s=1;  
    %%%%% uniformly distribute the injection energy with an energy bin
    for ii_e=1:length(mInj)       
        ind_e=ind_s+floor(Ne_s(mInj(ii_e)))-1;    % the end index for the energy bin ii_e
        Eps_in(ind_s:ind_e)=Emesh_injs(mInj(ii_e))+(rand(1,ind_e-ind_s+1)-0.5)*dE;
        ind_s=ind_e+1;                            % the start index for the energy bin ii_e+1
    end      
    Ne_s=Ne_s-floor(Ne_s);                        % update Ne_s, inject length(Eps_in) particles
    %%% input: the energy vector of the injected carriers, 
    %%% initialize the source injected carriers
    if Ntrs>0      
        xp_inj(1:Ntrs)=infs_x;   
        Ektp_inj(1:Ntrs)=max(Eps_in-Etops,1e-3); % initialize position and energy for particles above the SB
        vp_inj(1:Ntrs)=vE(Ektp_inj(1:Ntrs),Egh1,vF);
    end
        
    %%%%%%%%%%%%%% the drain injection
    Ne_d=Ne_d+Ninjd/sum(1./(1+exp((Emesh_injd-mu_d)./kBT)))*(1./(1+exp((Emesh_injd-mu_d)./kBT)));   % the source injection spectrum
    mInj=find(Ne_d>=1);  
    Ntrd=sum(floor(Ne_d(mInj)));  
    Ntrt=Ntrs+Ntrd; % total no. of transmitted particles
    ind_s=1;    
    for ii_e=1:length(mInj)       
        ind_e=ind_s+floor(Ne_d(mInj(ii_e)))-1;    % the end index for the energy bin ii_e
        Epd_in(ind_s:ind_e)=Emesh_injd(mInj(ii_e))+(rand(1,ind_e-ind_s+1)-0.5)*dE;
        ind_s=ind_e+1;      % the start index for the energy bin ii_e
    end      
    Ne_d=Ne_d-floor(Ne_d);
    %%% input: the energy vector of the injected carriers, 
    %%% initialize the drain injected carriers
    if Ntrd>0
        xp_inj((1:Ntrd)+Ntrs)=Lch-infs_x;   
        Ektp_inj((1:Ntrd)+Ntrs)=max(Epd_in-Etopd,1e-3); % initialize position and energy for particles above the SB    
        vp_inj((1:Ntrd)+Ntrs)=-vE(Ektp_inj((1:Ntrd)+Ntrs),Egh1,vF);   
    end   
    Ntrt=Ntrs+Ntrd;
    
    %%%%% add the injected particles to the 1D array data structure
    if Ntrt~=0
     %% add the injected carriers
        ind=length(xp)+1:length(xp)+Ntrt;
        xp(ind)=xp_inj;
        Ektp(ind)=Ektp_inj;
        vp(ind)=vp_inj;
        
    end
    %%%%%%%%%%% end of carrier injection %%%%%%%  