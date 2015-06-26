% Jing Guo, UFL

function [Ektpv xpv vpv Nh_bias,Idv, jSv, jDv, Emeshv]=mc(XI,E1,Vd_bias);
%%%%%% Numerical specification, inpute Em in the node grid, but output Ne,
%%%%%% Fn on the element grid
global q kBT 
global vF Egh1
% global set by inp_mc.m
global n t_step qsup Ntran Nstdy Degflag Balflag    % MC simulation parameters 
global Nkf Nxf                                      % E, k, and x grid
global Evec            % scattering rate
global  kfd    % the output using global parameters

dwire=n*1.42e-10*sqrt(3)/pi;
vF=9.8e5;               % m/s, Fermi velocity
flags.OutSflag=1e10;                                      % the velocity flag of the carriers that exit the source
flags.OutDflag=-1e10;                                     % the velocity flag of the carriers that exit the drain 
flags.IRflag=2;                                           % the velocity flag of the carriers that recombines 
flags.Degflag=Degflag;
flags.Balflag=Balflag;

%%%% conduction band Ec1 parameters
Ec1=E1;   % The conduction band edge.
mu_sc=0;                                             % the source Fermi level
mu_dc=-Vd_bias;                                      % the drain Fermi level
%%%%%% creat the energy mesh for carrier injection
%%% For Ec&Ev calculations, map valence band by flipping upside down
Eg_hyper=2*Egh1;  % hyperthetical bandgap
%%% the valence band holes treated by flipping the energy
Ev1=-(E1-Eg_hyper);
mu_sv=-mu_sc;                                             % the source Fermi level
mu_dv=-mu_dc;                                      % the drain Fermi level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ek=inline('Egh*(sqrt(1+(k./k0).^2)-1)','k','Egh','k0');
% vE=inline('vF*sqrt(1-(1./(1+E./Egh)).^2)','E','Egh','vF');
kE=inline('k0*sqrt(((E./Egh)+1).^2-1)','E','Egh','k0');   
kv=inline('k0*v./sqrt(vF^2-v.^2)','v','vF','k0');   

%%%%%%%%%%%%%% initialization %%%%%%%%%%%%%%%%%%
k0=2/(3*dwire);
kb=kE(max(Evec),Egh1,k0);                               % the boundary of k, [-kb, +kb]

Nnode=length(E1);
Fn_bias=zeros(Nnode,1);
Ne_bias=zeros(Nnode,1);
Lch=max(XI);                                        % the total channel length

%%%%% grid for collecting the distribution function f %%%%%%%
dxf=Lch/Nxf;                                        % the size of the position bin
xfd=(dxf/2):dxf:(Lch-dxf/2);                        % the position vector of the center of the bin
dkf=2*kb/Nkf;                                       % the size of the k bin

kfd=linspace(-kb+dkf/2,kb-dkf/2,Nkf);               % kfd normalized by k0
Npbin=2*2*(1/2/pi)*dxf*dkf/qsup;                    % the number of superparticle states in each x-k bin
dE=10e-3;

Eminc=min(Ec1);
Emaxc=max(Ec1)+20*kBT;
Emeshc=(Eminc+dE/2):dE:Emaxc; NEmeshc=length(Emeshc);    % energy grid for J(E)
%%%%%%% initialization %%%%%%%%%%%%%
%% injection carrier parameters
Ne_sc=[];                    % initialize the injection carrier spectrum
Ne_dc=[];
InOut_spec_c=zeros(NEmeshc,4);                              % the current spectrum                                          % no carrier at the beginning
fdc=zeros(Nxf,Nkf);                                  % initial distribution
sigmaNc=zeros(Nnode,1);                                % the number of particles summation
Ektpc=[]; xpc=[]; vpc=[];      % initialize the output parameters


%%%%%% creat the energy mesh for carrier injection
Eminv=min(Ev1);
Emaxv=max(Ev1)+20*kBT;
Emeshv=(Eminv+dE/2):dE:Emaxv; NEmeshv=length(Emeshv);    % energy grid for J(E)
%%%%%%% initialization %%%%%%%%%%%%%
%% injection carrier parameters
Ne_sv=[];                    % initialize the injection carrier spectrum
Ne_dv=[];
InOut_spec_v=zeros(NEmeshv,4);                              % the current spectrum                                          % no carrier at the beginning
fdv=zeros(Nxf,Nkf);                                  % initial distribution
sigmaNv=zeros(Nnode,1);                                % the number of particles summation
Ektpv=[]; xpv=[]; vpv=[];      % initialize the output parameters
 
%%%%%%%%%%% end of initialization %%%%%%%%%%%%%

%%%%%%%%%%% the Monte Carlo simulation %%%%%%%%%%%%%%
tic  
for ii_t=1:(Ntran+Nstdy)
    %%%%%%%%%%%%%%%% conduction band Ec1
    %%%% treat free flight and scattering
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%
    [xpc vpc Ektpc xpv vpv Ektpv]=scattering(XI, Ec1, flags, fdc, dwire,...
        xpc, vpc, Ektpc,xpv, vpv, Ektpv);  
    %% treat carriers out to contacts
    [Eps_out Epd_out xpc vpc Ektpc]=outflow(Ec1,flags,xpc,vpc,Ektpc);    
    %%% treat carrier injection
    [Eps_in Epd_in Ne_sc Ne_dc xpc vpc Ektpc]=inflow(XI,Ec1,Emeshc,mu_sc,mu_dc,Ne_sc,Ne_dc,xpc,vpc,Ektpc);    
    %%%% collect statistics for f(x,k), Ne, and current    
    if ii_t>Ntran  % steady state f(x,k), Ne, and I       
        [sigmaNc InOut_spec_c]=statQI(XI,xpc,sigmaNc,InOut_spec_c,Emeshc,Eps_in,Epd_in,Eps_out,Epd_out);
        if Degflag==1
            kp=kv(vpc,vF,k0);
            [fdt]=statfdt(xpc,kp,xfd,kfd,Npbin);
            fdc=fdc*(1-1/(ii_t-Ntran))+fdt*(1/(ii_t-Ntran));      % time average distribtion function
        end
    end  % end of % steady state f(x,k), Ne, and I
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% valance band Ev1
    %%%% treat free flight and scattering
    [xpv vpv Ektpv xpc vpc Ektpc]=scattering(XI, Ev1, flags, fdv, dwire,...
        xpv, vpv, Ektpv,xpc, vpc, Ektpc);
    %% treat carriers out to contacts
    [Eps_out Epd_out xpv vpv Ektpv]=outflow(Ev1,flags,xpv,vpv,Ektpv);    
    %%% treat carrier injection
    [Eps_in Epd_in Ne_sv Ne_dv xpv vpv Ektpv]=inflow(XI,Ev1,Emeshv,mu_sv,mu_dv,Ne_sv,Ne_dv,xpv,vpv,Ektpv);    
    %%%% collect statistics for f(x,k), Ne, and current    
    if ii_t>Ntran  % steady state f(x,k), Ne, and I       
        [sigmaNv InOut_spec_v]=statQI(XI,xpv,sigmaNv,InOut_spec_v,Emeshv,Eps_in,Epd_in,Eps_out,Epd_out);
        if Degflag==1
            kp=kv(vpv,vF,k0);
            [fdt]=statfdt(xpv,kp,xfd,kfd,Npbin);
            fdv=fdv*(1-1/(ii_t-Ntran))+fdt*(1/(ii_t-Ntran));      % time average distribtion function
        end
    end  % end of % steady state f(x,k), Ne, and I
end % end of ii_t
toc

%%% collect time average statistics  
%%% the charge density
u=zeros(Nnode,1); h=diff(XI);
u(2:(Nnode-1))=0.5*(h(1:(Nnode-2))+h(2:(Nnode-1))); u(1)=h(1)/2; u(Nnode)=h(Nnode-1)/2; %the element size

%%% output for Ec1
Ne_bias=(qsup/Nstdy)*(sigmaNc./u);     %Ne_bias=smooth((qsup/Nstdy)*(sigmaN./u),7);  
%%% the current spectrum 
jSc=(InOut_spec_c(:,1)-InOut_spec_c(:,3))*(1/Nstdy)*(q*qsup/t_step/dE);
jDc=(InOut_spec_c(:,4)-InOut_spec_c(:,2))*(1/Nstdy)*(q*qsup/t_step/dE);  
Idsc=sum(jSc)*dE
Iddc=sum(jDc)*dE
Idc=(Idsc+Iddc)/2;

%%% output for Ev1
Nh_bias=(qsup/Nstdy)*(sigmaNc./u);     %Ne_bias=smooth((qsup/Nstdy)*(sigmaN./u),7);  
%%% the current spectrum 
jSv=(InOut_spec_v(:,1)-InOut_spec_v(:,3))*(1/Nstdy)*(q*qsup/t_step/dE);
jDv=(InOut_spec_v(:,4)-InOut_spec_v(:,2))*(1/Nstdy)*(q*qsup/t_step/dE);
Idsv=sum(jSv)*dE
Iddv=sum(jDv)*dE
Idv=(Idsv+Iddv)/2;




