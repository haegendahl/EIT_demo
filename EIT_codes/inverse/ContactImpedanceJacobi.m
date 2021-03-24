
function Jz=ContactImpedanceJacobi(g,H,zimp,z,Eind,B)

%  - Compute Jacobian wrt contact impedances in EIT
%
% Synopsis: z is 1x3 cell containing electrode voltages z{1}, U.Current
% in z{2} and U.MeasField in z{3}, zimp contains contact impedances, Eind
% contains electrodes as rows wrt element indices, B contains boundary information
%
% Description:
%
% Examples:
%
% See also:
%

% - <optossav@physics.uku.fi>  -  16.06.2006
% Department of Physics, University of Kuopio
%
% Last modified - 20.06.2006 - (optossav)
% Modified from FemMatrix.m in EIDORS package

%% Compute Jacobian wrt contact impedances
Nel=max(size(zimp)); 
gN=max(size(g));
HN=max(size(H));
Jz=[];
[II1,C]=Current(Nel,gN,'adj');
C=C(:,1:Nel-1); %muutettu AN
%C = [ones(1,Nel-1);-eye(Nel-1)];
for jj=1:Nel
    %dA_dz=sparse(gN+Nel-1,gN+Nel-1);
    Dl=sparse(gN,Nel);
    Cl=sparse(gN,gN);
    s=zeros(Nel,1);
    for ii=1:length(Eind{jj}) %Would be better if Eind was a cell array
        ind=H(Eind{jj}(ii),:); % elektrodin jj elementin ii nodejen indeksit
        InE=jj;
        %bind=B{Eind{jj}(ii)};
        indAN1=Eind{jj};
        indAN2=indAN1(ii);
        bind=B(indAN2,:); % muutettu AN
        abc=g(bind,:);   
        bb2=triang1(abc);Bb2=zeros(max(size(ind)));

    eind=[find(bind(1)==ind),find(bind(2)==ind),find(bind(3)==ind)];
        Bb2(eind,eind)=bb2;
        Cl(ind,ind)=Cl(ind,ind)-1/(zimp(InE))^2*Bb2;
        s(InE)=s(InE)-1/(zimp(InE))^2*elektro(abc); 
        bb1=triang2(abc);Bb1=zeros(max(size(ind)),1);
        Bb1(eind)=bb1;    
        Dl(ind,InE)=Dl(ind,InE)+1/(zimp(InE))^2*Bb1;
    end %for ii
    El=diag(s);%sparse(diag(s));
    Dl=Dl*C;
    El=C'*El*C;
    dA_dz=[Cl,Dl;Dl.',El];
    DZ=-z{3}'*dA_dz*z{2};
    Jz=[Jz,DZ(:)];
end %for jj
  
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
