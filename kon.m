function konv=kon(D3d,D1d,knson,knsoff,rint,lDNA,varargin)

if lDNA==0
  if length(varargin)==0
    konv=knson/knsoff*pi./integral(@(q) 1./(D1d*q.^2+knsoff./(1+knson./(2*pi*D3d.*xK1divK0(abs(q)*rint)))),0,Inf);
  else
    kcap=varargin{1};
    konv=knson/knsoff*pi./integral(@(q) 1./(D1d*q.^2+knsoff./(1+knson./(2*pi*D3d.*xK1divK0(sqrt(q.^2+kcap/D3d)*rint)))),0,Inf);
  end
else
  if length(varargin)>0
    kcap=varargin{1};
  else
    kcap=findkcap(D3d,knson,rint,lDNA);
  end
  Wbulkcylq0=@(q) 1./(1+2*pi*D3d*xK1divK0(abs(q)*rint)/knson);
  Papps=@(s,a) a*exp(-a^2./(4*s))./(2*sqrt(pi*s.^3));
  lambdacoilq1=@(q) integral(@(s) Papps(s,4*sqrt(pi*D3d)*rint*lDNA).*Wbulkcylq0(sqrt(q^2+(s+kcap)/D3d)),0,Inf);
  konv=knson/knsoff*pi./integral(@(q) 1./(D1d*q.^2+knsoff*(1-arrayfun(lambdacoilq1,q))),0,Inf);
end
end

