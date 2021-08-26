% Semianalytical tool for bunch lengthening study 
% Bunch lengthing computation 
% author : Tianlong
% time   : 2021/08/26
clc;clear;
%% main parameters
% HEPS
% sigma_t = 40e-12;   % initial set rms bunch length  in ps
% sigma_E = 1.02e-3;  % natural energy spread
% alpha_c = 1.83e-5;  % momentum compaction factor
% C      = 1360.4;    % circumference in m
% h      = 756;       % harmonic number
% I0     = 200e-3;    % beam current in A
% U0     = 2.64e6;    % energy loss per turn in eV
% E0     = 6.0e9;     % beam energy in eV
% V_mc   = 3.176e6;   % MC voltage in V
% n_hc   = 3;         % harmonic number of HHC
% Q_hc   = 3.2e4;     % quality factor of HHC 
% R_hc   = Q_hc*47.65;% shunt impedance of HHC 
% fre_shift = 4300;   % detuning of HHC
% HALF
sigma_t = 10e-12;   % initial set rms bunch length  in ps
sigma_E = 6.45e-4;  % natural energy spread
alpha_c = 8.1e-5;  % momentum compaction factor
C      = 480;    % circumference in m
h      = 800;       % harmonic number
I0     = 300e-3;    % beam current in A
U0     = 198.8e3;    % energy loss per turn in eV
E0     = 2.2e9;     % beam energy in eV
V_mc   = 0.85e6;   % MC voltage in V
n_hc   = 3;         % harmonic number of HHC
Q_hc   = 5e5;     % quality factor of HHC 
R_hc   = Q_hc*90;% shunt impedance of HHC 
% fre_shift = 4300;   % detuning of HHC
fre_shift=detune_HC_calc(I0,n_hc,C,h,U0,V_mc,R_hc,Q_hc); % near-optimum lengthening condition
% pattern
%HEPS
% pattern = ones(1,h);
% pattern = zeros(1,h);pattern(1:680)=1;  % LONG GAP

% pattern = zeros(1,h);
% pat = [1:35,40:75,80:115,120:155,160:195,200:234,239:274,...
%     279:314,319:354,359:393,398:433,438:473,478:513,518:553,...
%     558:592,597:632,637:672,677:712,717:752]; % Distributed GAP
% pattern(pat)=1;

% pattern = zeros(1,h);
% pattern(1:170)=1; 
% pattern(190:359)=1;
% pattern(379:548)=1;
% pattern(568:737)=1;

% pattern = zeros(1,h);pattern(1:340)=1;pattern(379:718)=1;  % 

% HALF
pattern = ones(1,h); % complete fill 
% pattern=zeros(1,h);pattern(1:720)=1; % long gap fill

fillrate = length(find(pattern==1))/h;

% charge_ratio([1,25,30,55,60,85,90,115,120,145,150,174,179,204,209,234,239,264,269,294,299,324])=3.0;
% plot(charge_ratio);
% charge configuration
% 2% random charge distribution
charge_ratio = ones(1,h).*pattern;
% charge_ratio = charge_ratio+randn(1,h)*0.02.*pattern;
% charge_ratio = charge_ratio/sum(charge_ratio)*Bun_num;
% load('chargeratio1.mat');
HALF = machine(C,I0,U0,E0,sigma_t,sigma_E,alpha_c,h,V_mc,n_hc,R_hc,Q_hc,fillrate,fre_shift);
%% Iterative criterion setting
err = 1e-12; iteration = 0;  % for density distribution iteration
accuracy=1e-2;den_error = 1; % for synchrotron phase deviation iteration

%% save file
savename = ['HALF_Uniform_Current',num2str(I0*1e3),'mA_fre_',num2str(fre_shift),'_Hz.mat'];
%% HALF.fais_mc_whc = HALF.fais_nat;
HALF.qc    = charge_ratio.*pattern * HALF.qc;           % real charge per bunch
HALF.V_b  = HALF.qc * HALF.w_r * HALF.R_hc / HALF.Q_hc; 
%% initial distribution, assumed that all bunch are the same
tau = HALF.fai/HALF.w_rf; delta_tau = tau(2)-tau(1);

%% wake data (intrabunch motion)  data in column
tau_q = (0:length(tau)-1)'*delta_tau;
Wake_inter = -HALF.w_r *  R_hc /Q_hc*exp(-tau_q*HALF.w_r/2/Q_hc) .*cos(tau_q*HALF.w_r);
Wake_inter(1) = Wake_inter(1)/2;
% Wake_inter(:)=0;
%% BBR wake
% fr = 30e9;Rs = 3e3;
% Wake_BBR = 2*pi*fr*Rs*exp(-tau_q*2*pi*fr/2).*(cos(sqrt(0.75)*tau_q*2*pi*fr)-...
%     sin(sqrt(0.75)*tau_q*2*pi*fr)/2/sqrt(0.75));
% Wake_BBR(1) = Wake_BBR(1)/2;
% 
% Wake_inter = Wake_inter - Wake_BBR;
% plot(tau_q,Wake_inter);

%% initial distribution
Bun_num = length(find(pattern==1));
norm_den_dist = 1/(sqrt(2*pi)*sigma_t)*exp(-(tau/sqrt(2)/sigma_t).^2);
norm_den_dist_new = zeros(length(tau),Bun_num);
norm_den_dist = repmat(norm_den_dist',1,Bun_num);
TAU = repmat(tau',1,Bun_num);

%%
rot_decay_coef = 1i - 1 / (2 * HALF.Q_hc);               % rotating and decaying term
tau_list = zeros(1,h);
tau_list(2:h) = (1:h-1)*HALF.Tb;tau_list(1)=HALF.Tb*h;
V_load_coef = exp(rot_decay_coef * HALF.w_r * tau_list);
V_load_equi_coef = rot_decay_coef * HALF.w_r * HALF.T0;  % decaying coefficient per turn
exp_ang_coef   = -rot_decay_coef * HALF.w_r;
exp_angle      = exp(exp_ang_coef * TAU);  
%%
tic;
delta_fais = zeros(1,Bun_num);
MainIterationLoopCode;
% save data
save(savename,'HALF','norm_den_dist', 'delta_fais','BunchRMSLength','delta_tau','TAU','charge_ratio','Bun_num');
toc;
%%
cspeed = 299792458;
figure(1);
% recth1 =2:4:26;recth2= 260:4:284;recth=[recth1,recth2];
recth =[1:20:680];
% recth =[1,250,500];
for i=1:length(recth)
    TAUi=TAU(:,recth(i))+delta_fais(recth(i))/(2*pi*HALF.f_rf);
    zi  = TAUi*cspeed;
%     plot(zi,norm_den_dist(:,recth(i))/cspeed);xlim([-6e-2,6e-2]); hold on;
plot(TAUi,norm_den_dist(:,recth(i))*delta_tau,'Linewidth',1.5);xlim([-500e-12,500e-12]); hold on;

end
xlabel('\tau [s]');ylabel('norm density [s^-1]');
% legend('1th','35th','70th');
figure(2)
subplot(1,2,1)
plot(round(BunchRMSLength*1e12*1e6)/1e6,'.');hold on;
xlabel('bunch number');ylabel('\sigma_{\tau} [ps]');
set(gca,'FontName','Times New Roman','FontSize',12);
subplot(1,2,2)
Q_mean = sum(TAU.*norm_den_dist*delta_tau);
plot(round((delta_fais/HALF.w_rf+Q_mean)*1e12*1e6)/1e6,'.','Linewidth',1.5);hold on;
xlabel('bunch number');ylabel('<\tau> [ps]');
set(gca,'FontName','Times New Roman','FontSize',12);
%%