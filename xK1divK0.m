function out = xK1divK0(x)

K0=besselk(0,x);
out=x.*besselk(1,x)./K0;
out(K0==0)=x(K0==0);

end
