% PI structure KP(s+KI/s)
PI.m    = 5120/64; % ȡƽ����� 5120 , if 5120/2  32 tap FIR
% 64 tap FIR 
PI.tap  = 16;                            % for PI_Control_copy64
PI.dIQ  = PI.m/PI.tap;   % һ��IQ���     % for PI_Control_copy64
PI.d    = 800;  % �ӳټ��,�����Ǵ���0������ 800

PI.piIndex =0;
PI.Integral=0; % initial value PI ����ֵ
PI.RL      = HALF.R_mc; % ��ǻ�����迹

PI.Ts = PI.m * HALF.Tb;
% PI.Ts = PI.dIQ * HALF.Tb;
% PI.KP = 1; % PI �������ı���ϵ�� 0.2 1 2 
% PI.KI = HALF.w_rf/(2*HALF.Q_mc)*PI.Ts*PI.KP/64;  % PI �������Ļ���ϵ�� 10 50 100

PI.KP = 1; % PI �������ı���ϵ�� 0.2 1 2 
PI.KI = 0.002;  % PI �������Ļ���ϵ�� 10 50 100
% ���Է���KIҪ����С���Ա��⼤������������

PI.Ig0 = HALF.Ig_mc_0;                  % ���������ʸ��
PI.DIg = 0;
PI.Ig  = PI.Ig0;                        % ��ʼ Ig=Ig0
PI.Vg0 = HALF.Vg_mc_0;
PI.DI  = 0;
PI.Ig_track=[];

PI.RoverQ = HALF.R_mc/HALF.Q_mc;  % ��ǻ��R/Q������PI����
PI.Vrf_mc_track = zeros(1,PI.m);  % ����m��buckets������ǻѹʸ��

% RF phase modulation ���� to test
% mf ���ƶ�
PI.RFModulation_mf = 0.002;
% wm ����Ƶ��
PI.RFModulation_wm = 0.0012*2/HALF.T0*2*pi; %0.0025  0.002505
% rf frequency
PI.RFModulation_wf = HALF.w_rf;
PI.RFModulation_Tb = HALF.Tb;
PI.Detun_time      = HALF.det_angle_mc/HALF.w_rf;
