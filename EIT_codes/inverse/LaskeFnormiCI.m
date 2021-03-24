function [F,U] = LaskeFnormiCI(sigma1,U0,sigma0,sigma,args)

%sigma1 = sigma + step.*suunta;
ng = size(args{8},2);
I = find(sigma1 <= 0);
sigma1(I) = 1e-5;

sigma2 = args{8}*sigma1(1:ng,:);
z = sigma1(ng+1:end,:);

Uref = ForwardSolution3dCI(args{1},args{2},args{3},sigma2,z,args{5},'real');
%Uref = ForwardSolution3d2ndLinSigma2(args{1},args{2},args{9},args{3},10e3,sigma2,z,args{5},'real');

U = Uref.Electrode(:);

a = args{6}*(U0-U);
b = args{7}*(sigma1-sigma0);
F = norm(a)^2 + norm(b)^2;
