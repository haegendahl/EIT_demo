function [B,C,D,E] = genmatDfast_real(g,H,Element,Nel); 
  
NH = size(H,1);
Ng = size(g,1);
  
Br = zeros(4*NH,4);
Bc = zeros(4*NH,4);
Bv = zeros(4*NH,4);   

C = sparse(Ng,Ng);
D = sparse(Ng,Nel);
s = zeros(Nel,1);

a=0.58541020;b=0.13819660;
ip=[b b b;a b b;b a b;b b a];

L=[-1 1 0 0;-1 0 1 0;-1 0 0 1];
ss = ones(4,1);

pos = 1;
for ii=1:NH
  ind = H(ii,:);
  gg = g(ind,:);
  
  Jt = L*gg;
  iJt = inv(Jt);
  G = iJt*L;
  int = G'*G*abs(det(Jt))/6;

  id = ind(:);
  Br(pos:pos+3,:) = [id id id id];
  Bc(pos:pos+3,:) = [ind;ind;ind;ind];
  Bv(pos:pos+3,:) = int;
  
  if ~isempty([Element(ii).Electrode]), 
                                       
    bind=[Element(ii).Electrode{2}]'; 
    abc=g(bind,:);                       
    InE=Element(ii).Electrode{1};       
    s(InE) = s(InE) + elektro([abc]);

    eind=[find(bind(1)==ind),find(bind(2)==ind),find(bind(3)==ind)];
    bb1=triang2([abc]);Bb1=zeros(4,1);
    bb2=triang1([abc]);Bb2=zeros(4);

    Bb1(eind)=bb1;
    Bb2(eind,eind)=bb2;
    D(ind,InE) = D(ind,InE) - Bb1;  % minor slowdown, could be
                                         % indexed as well
    C(ind,ind) = C(ind,ind) + Bb2;
  end
  pos = pos + 4;
end
B = sparse(Br,Bc,Bv,Ng,Ng);

[II1 Chat] = Current(Nel,Ng,'adj');
Chat = sparse(Chat(:,1:Nel-1));

D = D*Chat;
E = sparse(diag(s));
E = Chat'*E*Chat;
