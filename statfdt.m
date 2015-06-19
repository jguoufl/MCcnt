%%% statistical function for the distribution function at a t step
function [fdt]=statfdt(xp,kp,xfd,kfd,Npbin)

Nxf=length(xfd);
dxf=xfd(2)-xfd(1);

%%%%%%%%% collect the distribution function
for ii_x=1:Nxf
    mx=find(abs(xp-xfd(ii_x))<dxf/2);
    fdt(ii_x,:)=(1/Npbin)*hist(kp(mx),kfd);
end