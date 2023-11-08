D3d=1e-1;
D1d=1e-4;
rint=3e-3;
knsoff=1e0;

lsl=sqrt(D1d/knsoff);

lsleff=@(knson) sqrt(knson/(2*pi*D3d))*lsl;
kon_react=@(knson) knson*2*lsl;
kon_diff=@(knson) 4*pi*D3d*lsleff(knson)./sqrt(log(lsleff(knson)/rint));

knsons=10.^(-2:0.1:4);
knson_reacts=knsons(1:21);
knson_diffs=knsons(11:end);

kons=[];
for i=1:length(knsons)
  kons(i)=kon(D3d,D1d,knsons(i),knsoff,rint,0);
end

kon_reacts=kon_react(knson_reacts);
kon_diffs=kon_diff(knson_diffs);

figure(2)
loglog(knsons,kons,'k-')
xlabel('{\it k}_{on}^{ns} ({\mu}m^2/ms)') %,'FontSize',14)
ylabel('{\it k}_{on} ({\mu}m^3/ms)') %,'FontSize',14)
ax=gca;
ax.FontSize=18;
refresh
hold on
  loglog(knson_reacts, kon_reacts,'r--')
  loglog(knson_diffs, kon_diffs,'b-.')
hold off
