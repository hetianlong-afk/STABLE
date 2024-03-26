% PI structure KP(s+KI/s)
PI.m    = 5120/64; % interval of buckets (number) to average the cavity voltage
PI.d    = 800;     % delay, in units of number of buckets

PI.piIndex =0;
PI.Integral=0; % initial value PI integral
PI.RL      = HALF.R_mc; % main cavity loaded shunt impedance

PI.Ts = PI.m * HALF.Tb;

PI.KP = 1;     % PI kp
PI.KI = 1e-4;  % PI ki

PI.Ig0 = HALF.Ig_mc_0;                  % generator current phasor
PI.DIg = 0;
PI.Ig  = PI.Ig0;                        % initial Ig=Ig0
PI.Vg0 = HALF.Vg_mc_0;
PI.DI  = 0;
PI.Ig_track=[];
PI.Ig_FBid =[];                         % used for PI feedback
PI.RoverQ = HALF.R_mc/HALF.Q_mc;  % main cavity R/Q
PI.Vrf_mc_track = zeros(1,PI.m);  % used for PI feedback

% RF phase modulation to do
% mf 调制度
PI.RFModulation_mf = 0.002;
% wm 调制频率
PI.RFModulation_wm = 0.0012*2/HALF.T0*2*pi; %0.0025  0.002505
% rf frequency
PI.RFModulation_wf = HALF.w_rf;
PI.RFModulation_Tb = HALF.Tb;
PI.Detun_time      = HALF.det_angle_mc/HALF.w_rf;
