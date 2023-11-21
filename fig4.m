D3d=1e-1;
rint=3e-3;
knson=1e4;
lDNA=1e3;

kcap=findkcap(D3d,knson,rint,lDNA)

Jcapu=@(u) knson*lDNA./u.^2./(1+knson./(2*pi*D3d*xK1divK0(sqrt(u/D3d)*rint)));
Jcap=@(t) talbot_inversion(Jcapu,t,8)';

ts=10.^(-16:.05:2);
Jcaps=[];
Regime1s=[];
Regime2and3s=[];
for i=1:length(ts)
  Jcaps(i)=Jcap(ts(i));
  Regime1s(i)=knson*lDNA*ts(i);
  Regime2and3s(i)=4*lDNA*rint*sqrt(pi*D3d*ts(i))+kcap*ts(i);
end

figure(4)
loglog(ts,Jcaps,'k-','LineWidth',2)
xlim([ts(1) ts(end)]);
xlabel('{\it t} (ms)')
ylabel('{\it J}_{cap}({\it t})')
ax=gca;
ax.FontSize=18;
refresh
hold on
  loglog(ts,Regime1s,'r--','LineWidth',2)
  loglog(ts,Regime2and3s,'b-.','LineWidth',2)
hold off

