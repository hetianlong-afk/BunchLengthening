function [delta_fais,X] = NewTon_calc(delta_fais,fais_mc_whc,U0,V_mc,V_load_comp_new,accuracy,HALF)
% 牛顿法迭代解同步相位
% 已知： delta_fais   U0  V_mc  V_load_comp_new
%  HALF.fais_mc_whc;

X_error = 1;
while X_error > accuracy
a = HALF.a;
coef1 = exp(a*(delta_fais'-delta_fais)).*V_load_comp_new;

sum_coef1_real = sum(real(coef1));
sum_coefa_real = sum(real(coef1'*a));

% matrix F 
F = -V_mc*sin(fais_mc_whc + delta_fais)+U0+sum_coef1_real;
% matrix A 
% A1= diag(V_mc*cos(HALF.fais_mc_whc + delta_fais));
% A2 = coef1*a;
% A3 = diag(a*sum_coef1_real);
% A = A1-A2-A3;
% X*A = -F
A = diag(V_mc*cos(fais_mc_whc + delta_fais))-real(coef1'*a)+diag(sum_coefa_real);

X = A\F';
delta_fais = delta_fais + X';
X_error = max(abs(X));
% figure(111)
% plot(delta_fais);
end