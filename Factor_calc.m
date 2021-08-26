function [F]=Factor_calc(sigma_t0,fre_hc,Q_hc,w_r)
tau=-5*sigma_t0:sigma_t0/1000:5*sigma_t0;
norm_den_dist=1/(sqrt(2*pi)*sigma_t0)*exp(-(tau/sqrt(2)/sigma_t0).^2);
fourier_tran=norm_den_dist.*exp(-1i*2*pi*fre_hc*tau+w_r*tau/2/Q_hc);
fourier_inte=(2*sum(fourier_tran)-fourier_tran(1)-fourier_tran(end))*(tau(2)-tau(1))/2;
F=abs(fourier_inte);
end