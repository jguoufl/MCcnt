function [sigmaN InOut_spec]=statQI(XI,xp,sigmaN,InOut_spec,Emesh,Eps_in,Epd_in,Eps_out,Epd_out)

%%% collect time average charge density, CIC approach
[ignore,k] = histc(xp,XI);        % k the index vecot XI(k)<=xp(mAll)<XI(k+1)
h=diff(XI);    % the spacing size vector of XI
pnr = (xp - XI(k))./h(k);     pnl = 1-pnr;       % CIC charge partition
for jj=1:length(xp)
    sigmaN(k(jj))=sigmaN(k(jj))+pnl(jj);
    sigmaN(k(jj)+1)=sigmaN(k(jj)+1)+pnr(jj);
end

%% current
InOut_spec(:,1)=InOut_spec(:,1)+hist(Eps_in,Emesh)';
InOut_spec(:,2)=InOut_spec(:,2)+hist(Epd_in,Emesh)';
InOut_spec(:,3)=InOut_spec(:,3)+hist(Eps_out,Emesh)';
InOut_spec(:,4)=InOut_spec(:,4)+hist(Epd_out,Emesh)';


