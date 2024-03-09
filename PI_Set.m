% PI structure KP(s+KI/s)
PI.m    = 5120/64; % 取平均间隔 5120 , if 5120/2  32 tap FIR
% 64 tap FIR 
PI.tap  = 16;                            % for PI_Control_copy64
PI.dIQ  = PI.m/PI.tap;   % 一组IQ间隔     % for PI_Control_copy64
PI.d    = 800;  % 延迟间隔,必须是大于0的整数 800

PI.piIndex =0;
PI.Integral=0; % initial value PI 积分值
PI.RL      = HALF.R_mc; % 主腔负载阻抗

PI.Ts = PI.m * HALF.Tb;
% PI.Ts = PI.dIQ * HALF.Tb;
% PI.KP = 1; % PI 控制器的比例系数 0.2 1 2 
% PI.KI = HALF.w_rf/(2*HALF.Q_mc)*PI.Ts*PI.KP/64;  % PI 控制器的积分系数 10 50 100

PI.KP = 1; % PI 控制器的比例系数 0.2 1 2 
PI.KI = 0.002;  % PI 控制器的积分系数 10 50 100
% 测试发现KI要尽量小，以避免激励束流纵向振荡

PI.Ig0 = HALF.Ig_mc_0;                  % 发射机电流矢量
PI.DIg = 0;
PI.Ig  = PI.Ig0;                        % 初始 Ig=Ig0
PI.Vg0 = HALF.Vg_mc_0;
PI.DI  = 0;
PI.Ig_track=[];

PI.RoverQ = HALF.R_mc/HALF.Q_mc;  % 主腔的R/Q，用于PI反馈
PI.Vrf_mc_track = zeros(1,PI.m);  % 保持m个buckets看到的腔压矢量

% RF phase modulation 参数 to test
% mf 调制度
PI.RFModulation_mf = 0.002;
% wm 调制频率
PI.RFModulation_wm = 0.0012*2/HALF.T0*2*pi; %0.0025  0.002505
% rf frequency
PI.RFModulation_wf = HALF.w_rf;
PI.RFModulation_Tb = HALF.Tb;
PI.Detun_time      = HALF.det_angle_mc/HALF.w_rf;
