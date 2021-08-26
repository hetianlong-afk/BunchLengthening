function [norm_den_dist,tau]=Norm_den_dist(Hf_wake,qc,intercoef,fai,V_mc,fais,n,k,faih,U0,frf,coef)
% fai :  phase
% V_mc:  Voltage of Main Cavity
% fais:  synthrotron phase
% k   :  Ratio of Harmonic Cavity Voltage
% U0  :  Energy loss per turn
% frf :  fundamental frequency
% coef:  potential coefficient
Hf = cos(fais)-cos(fais+fai)+k/n*(cos(faih)-cos(n*fai+faih))-U0*fai/V_mc;

% wake_conv = conv(normdendist',wake);
% Hf_wake = cumsum(wake_conv(1:length(wake)))*(deltatau^2);

den_dist = exp(coef*Hf+qc*intercoef*Hf_wake');

tau = fai/(2*pi*frf);
norm_den_dist = den_dist*2/(2*sum(den_dist)-den_dist(1)-...
    den_dist(end))/(tau(2)-tau(1));
end