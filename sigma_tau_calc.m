function [sigma_tau,mean_tau]=sigma_tau_calc(tau,norm_den_dist,flag)
% 计算束团RMS长度
% flag =1 表示输出数据；否则不输出
den_sum = sum(norm_den_dist);

mean_tau = sum(tau.*norm_den_dist)/den_sum;

sigma_tau = sqrt(sum(norm_den_dist.*(tau-mean_tau).^2)/den_sum);

if flag ==1
    disp(['RMS bunch length is ',num2str(sigma_tau*1e12),' [ps]']);
end
end