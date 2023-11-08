function kcap=findkcap(D3d,knson,rint,lDNA)

Wbulkcyl0u=@(u) 1./(1+2*pi*D3d*xK1divK0(sqrt(u/D3d)*rint)/knson);
Wbulkcyl=@(t) euler_inversion(Wbulkcyl0u,t,8)';
Jcapu=@(u) knson*lDNA./u.^2./(1+knson./(2*pi*D3d*xK1divK0(sqrt(u/D3d)*rint)));
Pforeignsurv=@(t) exp(-talbot_inversion(Jcapu,t,9)');
Papp=@(t,kcap) exp(-4*lDNA*rint*sqrt(pi*D3d*t)-kcap*t);
kcapeq=@(log2kcap) integral(@(t) Wbulkcyl(t).*(Pforeignsurv(t)-Papp(t,2^log2kcap)),0,Inf);
kcap=2^fzero(kcapeq,log2(D3d*lDNA));

end
