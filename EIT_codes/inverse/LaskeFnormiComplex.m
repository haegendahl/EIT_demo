function [F,U] = LaskeFnormiComplex(sigma1,U0,sigma0,sigma,args,Nel)

ns = length(sigma1)/2;
%sigma1 = sigma + step.*suunta;
%I = find(sigma1 <= 0);
%sigma1(I) = 1e-5;

sigma2 = reshape(args{8}*sigma1,ns,2);
sigma2 = complex(sigma2(:,1),sigma2(:,2));

Uref = ForwardSolution3d(args{1},args{2},args{3},sigma2,args{4},args{5},'comp');

%Ure = Uref.Electrode(1:Nel,:);
%Uim = Uref.Electrode(Nel+1:2*Nel,:);
U = Uref.Electrode(:); %[Ure(:);Uim(:)];

a = args{6}*(U0-U);
b = args{7}*(sigma1-sigma0);
F = norm(a)^2 + norm(b)^2;
