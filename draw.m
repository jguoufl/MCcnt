figure(1)   % current spectrum
plot(jS,Emesh,'linewidth',[2]); hold on
plot(jD,Emesh,'r--','linewidth',[2]); hold on
set(gca, 'fontsize',[24], 'linewidth',[2]);
set(gca,'position',[0.20 0.2 0.70 0.7]);
ylabel('E (eV)','fontsize',[24]);
xlabel('Je(A/eV)','fontsize',[24]);
%xlim([0 2e-4]);

% snapshot of carrier distribution at the end of the simulation
%load particles
mAll=find(abs(vp)>10);
figure(2)   % carrier distribution at the final bias point
plot(1e9*xp(mAll),Ektp(mAll)+interp1(XI,E1,xp(mAll)),'r.');   hold on
plot(1e9*XI,E1,'linewidth',[2]); hold on
xlim([0 max(XI)*1e9]);
set(gca, 'fontsize',[24], 'linewidth',[2]);
set(gca,'position',[0.15 0.2 0.74 0.7]);
xlabel('x (nm)','fontsize',[24]);
ylabel('E (eV)','fontsize',[24]);

% scattering rate vs. energy
figure(3)
plot(Evec,ga_tot,'linewidth',[2]); hold on
xlim([0 1]);
set(gca, 'fontsize',[24], 'linewidth',[2]);
set(gca,'position',[0.15 0.2 0.74 0.7]);
xlabel('E (eV)','fontsize',[24]);
ylabel('\Gamma (s^{-1})','fontsize',[24]);

