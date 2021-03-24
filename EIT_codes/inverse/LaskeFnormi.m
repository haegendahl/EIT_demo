function [F,U] = LaskeFnormi(sigma1,U0,sigma0,sigma,args)

%sigma1 = sigma + step.*suunta;
I = find(sigma1 <= 0);
sigma1(I) = 1e-5;

sigma2 = args{8}*sigma1;

Uref = ForwardSolution3d(args{1},args{2},args{3},1,sigma2,args{4},args{5},'real');
U = Uref.Electrode(:);

a = args{6}*(U0-U);
b = args{7}*(sigma1-sigma0);
F = norm(a)^2 + norm(b)^2;
