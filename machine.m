function [HALF]=machine(C,I0,U0,E0,tau_s,tau_z,sigma_t0,sigma_e0,alpha_c,h,V_mc,n_hc,R_hc,Q_hc,fillrate,fre_shift_hc,Q_mc,R_mc,fre_shift_mc)
% 创建机器参数，保存于HALF结构体中
cspeed = 299792458;
HALF.C       = C;
HALF.R       = C/(2*pi);
HALF.I0      = I0;
HALF.h       = h;
HALF.U0      = U0;
HALF.E0      = E0;
HALF.sigma_t0= sigma_t0;
HALF.sigma_e0= sigma_e0;
HALF.alpha_c = alpha_c;
HALF.tau_s   = tau_s;
HALF.tau_z   = tau_z;
HALF.V_mc    = V_mc;
HALF.n_hc    = n_hc;
HALF.R_hc    = R_hc;  
HALF.Q_hc    = Q_hc;
HALF.fillrate  = fillrate;
HALF.fre_shift_hc = fre_shift_hc;
HALF.R_mc    = R_mc;  
HALF.Q_mc    = Q_mc;
HALF.fre_shift_mc = fre_shift_mc;

HALF.gamma = E0/0.51099906/1e6;       % 电子静止能量  0.51099906 MeV
HALF.beta  = sqrt(1-1/(HALF.gamma^2));

HALF.T0 = HALF.C/cspeed/HALF.beta;    % 回旋周期
HALF.Tb = HALF.T0/HALF.h;             % 相邻 norminal bucket position 距离
HALF.f_rf = HALF.h/HALF.T0;           % 高频频率
HALF.w_rf = HALF.f_rf*2*pi;           % 高频角频率

HALF.fre_hc = HALF.f_rf * HALF.n_hc;  % 谐波腔未失谐谐振频率
HALF.wre_hc = HALF.fre_hc * 2*pi;     % 谐波腔未失谐谐振角频率
HALF.qc = HALF.T0*HALF.I0/HALF.h/HALF.fillrate;         % bunch 电荷量
HALF.wr_hc = HALF.w_rf*HALF.n_hc+HALF.fre_shift_hc*pi*2;% 谐波腔失谐谐振角频率
HALF.angle_hc = HALF.wr_hc * HALF.Tb;     % 谐波腔相邻 norminal bucket center 旋转角

HALF.fre_mc = HALF.f_rf;                  % 主腔未失谐谐振频率
HALF.wre_mc = HALF.fre_mc * 2*pi;         % 主腔未失谐谐振角频率
HALF.wr_mc = HALF.wre_mc + HALF.fre_shift_mc*pi*2;
HALF.angle_mc = HALF.wr_mc * HALF.Tb;     % 主腔相邻 norminal bucket center 旋转角

HALF.lambda_rads = HALF.T0 / HALF.tau_s;  % 辐射阻尼率   单位：圈数
HALF.lambda_radz = HALF.T0 / HALF.tau_z;  % 量子激发率   单位：圈数

HALF.fais_nat = pi-asin(HALF.U0/HALF.V_mc);% natural synchrotron phase
HALF.det_angle_hc= pi - atan(HALF.Q_hc*(HALF.wr_hc/HALF.wre_hc-HALF.wre_hc/HALF.wr_hc));
disp(['考虑谐波腔负载并假定均匀填充情况的失谐角度：',num2str(HALF.det_angle_hc/pi*180),' deg']);
HALF.fais_mc_whc = pi - asin((HALF.U0+2*HALF.I0*HALF.R_hc*cos(HALF.det_angle_hc)^2)/HALF.V_mc);
disp(['考虑谐波腔负载并假定均匀填充情况的同步相位：',num2str(HALF.fais_mc_whc),' rad  ',num2str(HALF.fais_mc_whc/pi*180),' deg']);
HALF.V_hc_load_0 = -2*HALF.I0*HALF.R_hc*cos(HALF.det_angle_hc)*exp(1i*(pi-HALF.det_angle_hc));
disp(['考虑谐波腔负载并假定均匀填充情况的初始负载：',num2str(HALF.V_hc_load_0)]);

HALF.det_angle_mc= atan(HALF.Q_mc*(HALF.wr_mc/HALF.wre_mc-HALF.wre_mc/HALF.wr_mc));
disp(['考虑主腔负载并假定均匀填充情况的失谐角度：',num2str(HALF.det_angle_mc/pi*180),' deg']);
HALF.V_mc_load_0 = 2*HALF.I0*HALF.R_mc*cos(HALF.det_angle_mc)*exp(1i*HALF.det_angle_mc);
disp(['考虑主腔负载并假定均匀填充情况的初始负载：',num2str(HALF.V_mc_load_0)]);

HALF.Vrf_ideal = -HALF.V_mc*exp(-1i*(pi/2-HALF.fais_mc_whc));
HALF.Vg_mc_init = HALF.Vrf_ideal-HALF.V_mc_load_0;
HALF.Vg_mc_track_0 = HALF.Vg_mc_init; % initial generator voltage phasor

HALF.Ig_mc_0 = HALF.Vg_mc_track_0/(2*HALF.R_mc*cos(HALF.det_angle_mc))*exp(-1i*HALF.det_angle_mc);
HALF.Vg_mc_0 = 2*pi*HALF.R_mc/HALF.Q_mc*HALF.Ig_mc_0;

HALF.drift_coef = HALF.alpha_c * HALF.T0 * (HALF.sigma_e0 / HALF.sigma_t0);
HALF.kick_coef  = 1 / (HALF.E0 * HALF.sigma_e0);

HALF.rfcoef1 = HALF.kick_coef * HALF.V_mc; % HALF.V_mc need be changed later

HALF.rfcoef2 = HALF.w_rf * HALF.sigma_t0;
% HALF.rfcoef2 = (HALF.w_rf-30e3*2*pi) * HALF.sigma_t0; % 改变主腔的谐振频率

HALF.radampcoef = 2 * HALF.lambda_rads;
HALF.quanexcoef = 2 * sqrt(HALF.lambda_radz);

HALF.ploss = HALF.kick_coef * HALF.U0;

HALF.VbImagFactor_hc   = 1/sqrt(4*Q_hc^2-1);
HALF.rot_coef_hc       = sqrt(1-1/(2*Q_hc)^2);

HALF.VbImagFactor_mc   = 1/sqrt(4*Q_mc^2-1);
HALF.rot_coef_mc       = sqrt(1-1/(2*Q_mc)^2);

% 计算负载角
HALF.angle_V_mc = Vb_angle_calc(real(HALF.Vrf_ideal),imag(HALF.Vrf_ideal));
HALF.angle_V_ge = Vb_angle_calc(real(HALF.Vg_mc_init),imag(HALF.Vg_mc_init));

HALF.angle_load = HALF.fais_mc_whc-(pi/2 - HALF.angle_V_ge-HALF.det_angle_mc);
disp(['初始负载角设置为：',num2str(HALF.angle_load/pi*180),'deg']);
end