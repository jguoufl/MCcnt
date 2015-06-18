function []=scattering(XI, E1, flags, fd, dwire)

global xp vp Ektp
global kBT q hbar
global t_step qsup Egh1 vF
global Evec ga_ap ga_op ga_tot tau_IR hw            % scattering rate
global kfd Nxf

OutSflag=flags.OutSflag;                                      % the velocity flag of the carriers that exit the source
OutDflag=flags.OutDflag;                                     % the velocity flag of the carriers that exit the drain 
IRflag=flags.IRflag;                                           % the velocity flag of the carriers that recombines 
Degflag=flags.Degflag;
Balflag=flags.Balflag;
infs=1e-3;   % in eV

Lch=max(XI); % the total channel length
dxf=Lch/Nxf;                                        % the size of the position bin

%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ek=inline('Egh*(sqrt(1+(k./k0).^2)-1)','k','Egh','k0');
vE=inline('vF*sqrt(1-(1./(1+E./Egh)).^2)','E','Egh','vF');
kE=inline('k0*sqrt(((E./Egh)+1).^2-1)','E','Egh','k0');   
kv=inline('k0*v./sqrt(vF^2-v.^2)','v','vF','k0'); 

k0=2/(3*dwire);
Egh1=0.42/(dwire*1e9);                              % the half band gap
kb=kE(max(Evec),Egh1,k0);                               % the boundary of k, [-kb, +kb]
dkf=kfd(2)-kfd(1);

%%%% simulation of free fly and scattering
Np_max=length(xp);
%% initialization
    noAll=Np_max;  % initialize the number of carriers to be treated for free fly& scattering
    mAll=1:noAll; % initizlize index of carriers for free fly & scattering
    t=zeros(1,Np_max);   % initialize time for the subensemble mAll       
    
    %% iteration to treat carriers till the end of the time step
    while noAll>0             
        %%%%%%%%% (1) determine the scattering time %%%%%%%%%%%%%%%%%%%%%%%%
        EktI=Ektp(mAll);                            % the initial kinetic energy
        E1I=interp1(XI,E1,xp(mAll),'nearest','extrap');                % the initial potential energy
        % propose the end-of-step position of all particles
        Xend=xp(mAll)+vp(mAll).*(t_step-t(mAll));              
        % compute the proposed end-of-step kinetic energy of all particles
        E1F=interp1(XI,E1,Xend,'nearest','extrap'); % the final potential energy
        EktF=max(EktI+E1I-E1F,zeros(1,noAll));      % the final kinetic energy
        ga_I=interp1(Evec,ga_tot,EktI);
        ga_F=interp1(Evec,ga_tot,EktF);
        ga=zeros(1,Np_max);                         % the cap of the total scattering rate
        ga(mAll)=0.5*(1+sign((EktI-hw).*(EktF-hw))).*max(ga_I,ga_F)...
            +0.5*(1-sign((EktI-hw).*(EktF-hw))).*max(ga_tot);
        %%%%%% update time at the end of scattering       
        if Balflag==1
            t_prime=t_step*ones(1,noAll);           % ballistic transport
        else
            t_prime=min(t_step*ones(1,noAll),t(mAll)+(1./ga(mAll)).*log(1./rand(1,noAll))); 
        end
        tInc=t_prime-t(mAll);                       % update free flight time;
        t(mAll)=t_prime;                            % update the clock
        %%%%%%%%%%%%%(2.1) treat the thermal current
        xfly=zeros(1,Np_max);
        xfly(mAll)=xp(mAll)+tInc.*vp(mAll);         % the position at the end of the free flight
                       
        mExitD=mAll(find(xfly(mAll)>=Lch));         % the index of the subensemble that exit the drain;
        vp(mExitD)=OutDflag;
        mExitS=mAll(find(xfly(mAll)<=0));           % the index of the subensemble that exit the source;
        vp(mExitS)=OutSflag;               
             
        mRest1=mAll(find((xfly(mAll)>0)&(xfly(mAll)<Lch)));  % the index of the particles in the device
        if length(mRest1)>0       
            %%%%%%% (2.2) treat reflection by barriers
            Ekt_mR1=(interp1(XI,E1,xp(mRest1))+Ektp(mRest1))-interp1(XI,E1,xfly(mRest1)); % the kinetic energies of mRest1
            mBedge=mRest1(find(Ekt_mR1<=0));  % the reflected particle index in the top ensemble
            vp(mBedge)=-vp(mBedge);                         % treat the reflected carriers
                      
            %%%%%(2.3) free flight and scattering %%%%
            %%% treat free flight
            mRest2=mRest1(find(Ekt_mR1>0));                 % free flight & scattering carriers
            Ektp(mRest2)=Ekt_mR1(Ekt_mR1>0);
            xp(mRest2)=xfly(mRest2);
            vp(mRest2)=sign(vp(mRest2)).*vE(Ektp(mRest2),Egh1,vF);
                    
            % find the scattered subensemble in mRest2, 
            mScat=mRest2(find(t(mRest2)<t_step));
            mAll=mScat;
            noAll=length(mAll);
              
            %%%% treat scattering  %%%%%%%%%
            if noAll>0

                %%%% compute scattering rates for all scattered particles
                ga_ap_mScat=interp1(Evec,ga_ap,Ektp(mScat)); % AP scattering rate
                ga_op_mScat=interp1(Evec,ga_op,Ektp(mScat)); % OP scattering rate
                ga_IR_mScat=1/tau_IR*ones(1,length(mScat));  % IR emission rate
                ga_tot_mScat=ga(mScat);                      % total scattering rate
                ga_self_mScat=ga_tot_mScat-ga_ap_mScat-ga_op_mScat; % the self scattering rate
 
                rand_sca=rand(1,noAll);
                mAp=mScat(find(rand_sca<=(ga_ap_mScat./ga_tot_mScat)));                 % AP scattering
                mOp=mScat(find((rand_sca>(ga_ap_mScat./ga_tot_mScat)) &...
                    (rand_sca<=((ga_ap_mScat+ga_op_mScat)./ga_tot_mScat))));            % OP scattering
                mIR=mScat(find((rand_sca>((ga_ap_mScat+ga_op_mScat)./ga_tot_mScat)) &...
                    (rand_sca<=((ga_ap_mScat+ga_op_mScat+ga_IR_mScat)./ga_tot_mScat)))); % IR emission
                                                 
                % 2.3.1 AP scattering 
                if (length (mAp)>0)   % AP scattering 
                    if Degflag==1   % Pauli exclusion
                        ki_ap=kv(vp(mAp),vF,k0);
                        kf_ap=-ki_ap;   
                        fd_ap=diag(fd(ceil(xp(mAp)*(1/dxf)),ceil((kf_ap+kb)*(1/dkf))));
                        mAp=mAp(find(rand(length(mAp),1)>fd_ap));       % rejection technique for degeneracy
                    end    
                     % full-band AP back scattering
                    if length(mAp)>0     
                        ki_ap=kv(vp(mAp),vF,k0);
                        kf_ap=-ki_ap;   
                        Ektp(mAp)=Ek(kf_ap,Egh1,k0);    % carrier scatters to just above E1
                        vp(mAp)=sign(kf_ap).*vE(Ektp(mAp),Egh1,vF);   % backscattering                       
                    end                 
                end                                         
                % 2.3.2 OP scattering 
                if (length (mOp)>0)   
                    if Degflag==1   % Pauli exclusion
                        Ekt_op=max(infs, Ektp(mOp)-hw);               % the projected final kinetic energy 
                        kp_op=-sign(vp(mOp)).*kE(Ekt_op,Egh1,k0);
                        xp_op=xp(mOp);
                        fd_op=diag(fd(ceil(xp_op*(1/dxf)),ceil((kp_op+kb)*(1/dkf))));
                        mOp=mOp(find(rand(length(mOp),1)>fd_op));       % rejection technique for degeneracy                    
                    end                                       
                    % OP back scattering
                    if length(mOp)>0                        
                        %Nop_x=Nop_x+hist(xp(mOp),XI)';  
                        Ektp(mOp)=max(infs, Ektp(mOp)-hw);            % carrier kinetic energy
                        vp(mOp)=-sign(vp(mOp)).*vE(Ektp(mOp),Egh1,vF);  % carrier velocity                          
                    end 
                end     % if (length (mOp)>0) 
                 % 2.3.3 IR emission in mIR
                 if (length (mIR)>0)   % IR emission
                    if Degflag==1   % apply rejection rule for treating degeneracy
                        kp_IR=kv(vp(mIR),vF,k0);
                        fd_IR=diag(fd(Nxf+1-ceil(xp(mIR)*(1/dxf)),ceil((kp_IR+kb)*(1/dkf))));
                        mIR=mIR(find(rand(length(mIR),1)<fd_IR));   % IR emission, in mAll
                    end               
                    % IR emission
                    vp(mIR)=IRflag; % anialate the recombined carriers
                    mAll=mAll(vp(mAll)~=IRflag);  % exclude the recombined carriers
                    noAll=length(mAll);
            
                    if ii_t>Ntran   % light spot position and spectrum statistics
                        Nx_IR=Nx_IR+hist(xp(mIR),xfd)';
                        NE_IR=NE_IR+hist(2*Ektp(mIR),Emesh)';
                    end
                end      % end of if (length (mIR)>0)   
            end         % if noAll>0
        end             % if length(mRest1)>0                          
    end                 % while noAll>0
    %%%%%%%%%%% end of treating carrier free fly and scattering
    