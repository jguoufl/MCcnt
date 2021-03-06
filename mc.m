% Jing Guo, UFL

function [Ektp xp vp Ne_bias,Id, jS, jD, Emesh]=mc(XI,E1,Vd_bias);
%%%%%% Numerical specification, inpute Em in the node grid, but output Ne,
%%%%%% Fn on the element grid
global q kBT 
global vF Egh1
% global set by inp_mc.m
global n t_step qsup Ntran Nstdy Degflag Balflag    % MC simulation parameters 
global Nkf Nxf                                      % E, k, and x grid
global Evec            % scattering rate
%%% output global
global  kfd fdn   % the output using global parameters

dwire=n*1.42e-10*sqrt(3)/pi;
vF=9.8e5;               % m/s, Fermi velocity
flags.OutSflag=1e10;                                      % the velocity flag of the carriers that exit the source
flags.OutDflag=-1e10;                                     % the velocity flag of the carriers that exit the drain 
flags.IRflag=2;                                           % the velocity flag of the carriers that recombines 
flags.Degflag=Degflag;
flags.Balflag=Balflag;

%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ek=inline('Egh*(sqrt(1+(k./k0).^2)-1)','k','Egh','k0');
% vE=inline('vF*sqrt(1-(1./(1+E./Egh)).^2)','E','Egh','vF');
kE=inline('k0*sqrt(((E./Egh)+1).^2-1)','E','Egh','k0');   
kv=inline('k0*v./sqrt(vF^2-v.^2)','v','vF','k0');   

%%%%%%%%%%%%%% initialization %%%%%%%%%%%%%%%%%%
Ektp=[]; xp=[]; vp=[];      % initialize the output parameters
k0=2/(3*dwire);
Egh1=0.42/(dwire*1e9);                              % the half band gap
kb=kE(max(Evec),Egh1,k0);                               % the boundary of k, [-kb, +kb]

Nnode=length(E1);
Fn_bias=zeros(Nnode,1);
Ne_bias=zeros(Nnode,1);
Lch=max(XI);                                        % the total channel length
mu_s=0;                                             % the source Fermi level
mu_d=-Vd_bias;                                      % the drain Fermi level

%%%%%% creat the energy mesh for carrier injection
dE=10e-3;
[Emin, Imin]=min(E1);
Emax=max(E1)+20*kBT;
Emesh=(Emin+dE/2):dE:Emax; NEmesh=length(Emesh);    % energy grid for J(E)
%% injection carrier parameters
Ne_s=[];                    % initialize the injection carrier spectrum
Ne_d=[];

%%%%% grid for collecting the distribution function f %%%%%%%
dxf=Lch/Nxf;                                        % the size of the position bin
xfd=(dxf/2):dxf:(Lch-dxf/2);                        % the position vector of the center of the bin
dkf=2*kb/Nkf;                                       % the size of the k bin

kfd=linspace(-kb+dkf/2,kb-dkf/2,Nkf);               % kfd normalized by k0
Npbin=2*2*(1/2/pi)*dxf*dkf/qsup;                    % the number of superparticle states in each x-k bin
%%%%%%% initialization %%%%%%%%%%%%%
InOut_spec=zeros(NEmesh,4);                              % the current spectrum
                                          % no carrier at the beginning
fd=zeros(Nxf,Nkf);                                  % initial distribution
sigmaN=zeros(Nnode,1);                                % the number of particles summation
Nx_IR=zeros(Nxf,1);                                 % the emission position
NE_IR=zeros(NEmesh,1);                              % the emission energy spectrum
 
%%%%%%%%%%% end of initialization %%%%%%%%%%%%%
%%%%%%%%%%% the Monte Carlo simulation %%%%%%%%%%%%%%
tic  
for ii_t=1:(Ntran+Nstdy)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% treat free flight and scattering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [xp vp Ektp]=scattering(XI, E1, flags, fd, dwire, xp, vp, Ektp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% treat carriers out to contacts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [Eps_out Epd_out xp vp Ektp]=outflow(E1,flags,xp,vp,Ektp);
    
    %%% treat carrier injection
    [Eps_in Epd_in Ne_s Ne_d xp vp Ektp]=inflow(XI,E1,Emesh,mu_s,mu_d,Ne_s,Ne_d,xp,vp,Ektp);
    
    %%%% collect statistics for f(x,k), Ne, and current    
    if ii_t>Ntran  % steady state f(x,k), Ne, and I
        
        [sigmaN InOut_spec]=statQI(XI,xp,sigmaN,InOut_spec,Emesh,Eps_in,Epd_in,Eps_out,Epd_out);
        if Degflag==1
            kp=kv(vp,vF,k0);
            [fdt]=statfdt(xp,kp,xfd,kfd,Npbin);
            fd=fd*(1-1/(ii_t-Ntran))+fdt*(1/(ii_t-Ntran));      % time average distribtion function
        end
    end  % end of % steady state f(x,k), Ne, and I
end % end of ii_t
toc

%%% collect time average statistics  
%%% the charge density
u=zeros(Nnode,1); h=diff(XI);
u(2:(Nnode-1))=0.5*(h(1:(Nnode-2))+h(2:(Nnode-1))); u(1)=h(1)/2; u(Nnode)=h(Nnode-1)/2; %the element size
Ne_bias=(qsup/Nstdy)*(sigmaN./u);     %Ne_bias=smooth((qsup/Nstdy)*(sigmaN./u),7);  
%%% the current spectrum
jS=(InOut_spec(:,1)-InOut_spec(:,3))*(1/Nstdy)*(q*qsup/t_step/dE);
jD=(InOut_spec(:,4)-InOut_spec(:,2))*(1/Nstdy)*(q*qsup/t_step/dE);
%Id=0.5*sum(jD+jS)*dE   
Ids=sum(jS)*dE
Idd=sum(jD)*dE
Id=(Ids+Idd)/2;




