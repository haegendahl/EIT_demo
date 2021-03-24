function int=tetra2ndLinSigma(g,sigma,w,ip,np);

int=0;
for ii=1:np
  k=ip(ii,1);k4=4*k;
  n=ip(ii,2);n4=4*n;
  c=ip(ii,3);c4=4*c;
  
  L=[k4-3+n4+c4 k4-1 0 0 -8*k+4-n4-c4 n4 -n4 -c4 c4 0;
        n4-3+k4+c4 0 n4-1 0 -k4 k4 -8*n+4-k4-c4 -c4 0 c4;
        c4-3+k4+n4 0 0 c4-1 -k4 0 -n4 -8*c+4-k4-n4 k4 n4];
  S=[1-k-n-c;k;n;c]'*sigma;
  Jt=L*g;
  dJt=abs(det(Jt));
  G=inv(Jt)*L;
  int=int+w(ii)*S*dJt*G'*G;
end
