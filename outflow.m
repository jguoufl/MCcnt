%%% outflow from the device to the contacts

function [Eps_out Epd_out xp vp Ektp]=outflow(E1,flags,xp,vp,Ektp)
OutSflag=flags.OutSflag;                                      % the velocity flag of the carriers that exit the source
OutDflag=flags.OutDflag;                                     % the velocity flag of the carriers that exit the drain 
IRflag=flags.IRflag; 

%% find energy distributions of carriers out to S and D
mExitS=find(vp==OutSflag); % out to source
Eps_out=Ektp(mExitS)+E1(1);
mExitD=find(vp==OutDflag);  % out to drain
Epd_out=Ektp(mExitD)+E1(end);
%% eliminate carriers out to contacts
Out_I=find((vp==OutSflag)|(vp==OutDflag)|(vp==IRflag)); %mark carriers exited
 xp(Out_I)=[]; Ektp(Out_I)=[]; vp(Out_I)=[];  % eliminate carriers
  