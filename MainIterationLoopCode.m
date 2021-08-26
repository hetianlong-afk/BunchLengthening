while den_error > err && iteration<=1000 % the total number of iterations
iteration = iteration + 1;
% 当前圈负载电压矢量
V_load_init = HALF.V_b(pattern==1) .* sum(norm_den_dist.*exp_angle) * delta_tau;
% N圈之后负载电压矢量
V_load_equi = V_load_init / (1-exp(V_load_equi_coef)); %%% notice !!!
V_load_comp = zeros(Bun_num,h);
jj=0;
for i=1:h
    if pattern(i)==1
        jj=jj+1;
        V_load_equi_ith = V_load_equi(jj)*V_load_coef;
        V_load_comp(jj,i:end) = V_load_equi_ith(1:end-i+1);
        if i>1
            V_load_comp(jj,1:i-1)= V_load_equi_ith(end-i+2:end);
        end
    end
end
V_load_comp_new = V_load_comp(:,pattern==1);
V_load_final  = sum(V_load_comp_new);
%% 迭代求解同步相位fais  accuracy=2e-3;
delta_fais = zeros(1,Bun_num);
[delta_fais,X] = NewTon_calc(delta_fais,HALF.fais_mc_whc,U0,V_mc,V_load_comp_new,accuracy,HALF);
%%
% delta_fais = zeros(1,Bun_num);
% V_load_sum = sum(V_load_comp.*exp(exp_ang_coef*delta_fais'/HALF.w_rf)); % exp_ang_coef 符号为正
% V_load_new = V_load_sum(pattern==1).*exp(-exp_ang_coef*delta_fais/HALF.w_rf);

V_load_new = sum(V_load_comp_new.*exp(HALF.a*(delta_fais'-delta_fais)));

% k factor and detuning angle
eof = 1e-50;
Re = real(V_load_new);
Im = imag(V_load_new)+eof;
[Vb_angle_new]=round(Vb_angle_calc(Re,Im)*1e16)/1e16;
k_factor_new = round(sqrt(Re.^2+Im.^2)/V_mc*1e16)/1e16;

% [Vb_angle_new]=Vb_angle_calc(Re,Im);
% k_factor_new = sqrt(Re.^2+Im.^2)/V_mc;

if mod(iteration,50) ==0
    figure(12);
    subplot(2,1,1);
%     title('Amplitude ratio');
    plot(k_factor_new);hold on;
    xlabel('Bunch number');
    ylabel('Amplitude ratio');
    subplot(2,1,2);
%     title('Phase [rad]');
    plot(Vb_angle_new);hold on;
    xlabel('Bunch number');
    ylabel('Phase [rad]');
end

% synchrotron phase
fais = delta_fais + HALF.fais_mc_whc;              % 新的同步相位
nfaih = pi/2 - Vb_angle_new;
BunchRMSLength=zeros(1,Bun_num);
% 束团内尾势
wake_conv = conv2(norm_den_dist,Wake_inter);
Hf_wake = cumsum(wake_conv(1:length(Wake_inter),:))*(delta_tau^2); % 每列表示一个束团分布

j=0;
for i = 1:h
    if pattern(i)==1
        j=j+1;
        [norm_den_dist_new(:,j),TAU(:,j)] = Norm_den_dist(Hf_wake(:,j),HALF.qc(i),...
            HALF.poten_inter_coef,HALF.fai,HALF.V_mc,fais(j),HALF.n_hc,...
            k_factor_new(j),nfaih(j),HALF.U0,HALF.f_rf,HALF.poten_coef);
        norm_den_dist_ith = norm_den_dist_new(:,j);
        tau_ith = TAU(:,j);
        [BunchRMSLength(j),~]=sigma_tau_calc(tau_ith',norm_den_dist_ith',0);
    end
end

den_error = max(max(abs(norm_den_dist_new-norm_den_dist)./(norm_den_dist+1e-20)));
r=rand(1)/2;norm_den_dist =(1-r).*norm_den_dist + r.* norm_den_dist_new;
% norm_den_dist =0.5*(norm_den_dist + norm_den_dist_new);
if mod(iteration,100)==0
    disp(['iteration = ',num2str(iteration)]);
    figure(22);
    plot(delta_fais);
    figure(23);
    plot(TAU(:,1),norm_den_dist(:,1));hold on;
    figure(24);
subplot(1,2,1)
plot(round(BunchRMSLength*1e12*1e6)/1e6,'.');hold on;
xlabel('bunch number');ylabel('\sigma_{\tau} [ps]');
set(gca,'FontName','Times New Roman','FontSize',12);
subplot(1,2,2)
Q_mean = sum(TAU.*norm_den_dist*delta_tau);
plot(round((delta_fais/HALF.w_rf+Q_mean)*1e12*1e8)/1e8,'.','Linewidth',1.5);hold on;
xlabel('bunch number');ylabel('<\tau> [ps]');
set(gca,'FontName','Times New Roman','FontSize',12);
end
end