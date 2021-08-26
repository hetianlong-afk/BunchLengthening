function [HALF]=machine(C,I0,U0,E0,sigma_t,sigma_E,alpha_c,h,V_mc,n_hc,R_hc,Q_hc,fillrate,fre_shift)
% machine parameter structure
cspeed = 299792458;   
HALF.C       = C;              % circumference
HALF.R       = C/(2*pi);       % 
HALF.I0      = I0;             % beam current
HALF.h       = h;              % harmonic number
HALF.U0      = U0;             % energy loss per turn
HALF.E0      = E0;             % beam energy
HALF.sigma_t0= sigma_t;        % initial rms length in ps
HALF.sigma_z0= sigma_t*cspeed; % initial rms length in m
HALF.sigma_E0= sigma_E;        % natrual energy spread
HALF.alpha_c = alpha_c;        % momentum compaction factor
HALF.V_mc    = V_mc;           % MC voltage
HALF.n_hc    = n_hc;           % harmonic number of HHC 
HALF.R_hc    = R_hc;           % Shunt impedance of HHC 
HALF.Q_hc    = Q_hc;           % quality factor of HHC 
HALF.fillrate  = fillrate;     % filling rate of buckets
HALF.fre_shift = fre_shift;    % detuning of HHC 
HALF.T0 = HALF.C/cspeed;       % revolution time
HALF.Tb = HALF.T0/HALF.h;      % reference phase interval of two adjacent buckets
HALF.f_rf = HALF.h/HALF.T0;    % main frequency
HALF.fre_hc = HALF.f_rf * HALF.n_hc; % harmonic frequency
HALF.w_rf = HALF.f_rf*2*pi;          % main angular frequency
HALF.wre_hc = HALF.fre_hc * 2*pi;    % harmonic angular frequency
HALF.qc = HALF.T0*HALF.I0/HALF.h/HALF.fillrate;    % average chrage per bunch
HALF.w_r = HALF.w_rf*HALF.n_hc+HALF.fre_shift*pi*2;% real angular frequency of HHC 

% HALF.n_hc=HALF.w_r/HALF.w_rf;% should not be a integral
HALF.a = -1i*HALF.w_r/HALF.w_rf;
HALF.Factor = Factor_calc(120e-12,HALF.fre_hc,HALF.Q_hc,HALF.w_r);

HALF.V_b = HALF.qc * HALF.w_r * HALF.R_hc / HALF.Q_hc;
HALF.angle = HALF.w_r * HALF.Tb;        

fai_max = HALF.h * HALF.sigma_z0 / HALF.R*20;
HALF.fai = linspace(-fai_max,fai_max,3e2+1); % interested phase range
HALF.fais_nat = pi-asin(HALF.U0/HALF.V_mc);  % natural synchrotron phase
HALF.det_angle= pi - atan(HALF.Q_hc*(HALF.w_r/HALF.wre_hc-HALF.wre_hc/HALF.w_r));
disp(['detuning angle in uniform complete fill:',num2str(HALF.det_angle/pi*180),' deg']);
HALF.fais_mc_whc = pi - asin((HALF.U0+2*HALF.I0*HALF.Factor*HALF.R_hc*cos(HALF.det_angle)^2)/HALF.V_mc);
disp(['synchrotron phase in uniform complete fill:',num2str(HALF.fais_mc_whc),' rad']);
HALF.V_load_0 = -2*HALF.I0*HALF.R_hc*cos(HALF.det_angle)*exp(1i*(pi-HALF.det_angle));
disp(['initial beam loading voltage phasor:',num2str(HALF.V_load_0)]);
HALF.poten_coef = HALF.V_mc/(2*pi * HALF.h * HALF.alpha_c * HALF.E0 * HALF.sigma_E0^2);% 势能系数
HALF.poten_inter_coef = 1/(HALF.E0 * HALF.T0 * HALF.alpha_c * HALF.sigma_E0^2); %束内尾势系数
end