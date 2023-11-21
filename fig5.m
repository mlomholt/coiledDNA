D3d=1e-1;
D1d=1e-4;
rint=3e-3;
knsoff=1e0;
knson=1e4;

lsl=sqrt(D1d/knsoff);
lsleff=sqrt(knson/(2*pi*D3d))*lsl;

a=@(lDNA) 4*sqrt(pi)*rint*lDNA*sqrt(D3d);
kon_diff_coil=@(lDNA) 4*pi*D3d*lsleff.*sqrt(1./log(1./sqrt(lsleff^(-2)+lDNA)/rint)+4*rint^2*lDNA*log(1+knson/(2*pi*D3d)));

lDNAs=10.^(-10:0.2:-3)*1e6;
kons=[];
konapps=[];
koncyls=[];
kon_diff_coils=[];
kcaps=[];
for i=1:length(lDNAs)
  kcap=findkcap(D3d,knson,rint,lDNAs(i));
  kcaps(i)=kcap;
  kons(i)=kon(D3d,D1d,knson,knsoff,rint,lDNAs(i),kcap);
  koncyls(i)=kon(D3d,D1d,knson,knsoff,rint,0);
  kon_diff_coils(i)=kon_diff_coil(lDNAs(i));
  [i length(lDNAs) lsleff*sqrt(lDNAs(i)) kons(i)/koncyls(i) kon_diff_coils(i)/kon_diff_coil(0) lDNAs(i)*D3d/kcap]
end

figure(5)
semilogx(sqrt(lDNAs)*lsleff,kons./koncyls,'k-','LineWidth',2)
xlabel('{\it l}_{DNA}^{1/2}{\it l}_{sl}^{eff}') 
ylabel('{\it k}_{on}/{\it k}_{on}^{\rm cyl}') 
yline(0)
xlim([sqrt(lDNAs(1))*lsleff sqrt(lDNAs(end))*lsleff]);
ax=gca;
ax.FontSize=18;
refresh
hold on
  semilogx(sqrt(lDNAs)*lsleff,kon_diff_coils/kon_diff_coil(0),'r--','LineWidth',2)
hold off

