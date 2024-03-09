clc;clear;
cspeed = 299792458;
Ib     = 260e-3;
h      = 800;
C      = 480;
R      = C/(2*pi);

U0     = 198.8e3 ;   % energy loss per turn 218e3
V_mc   = 0.85e6; % 主腔腔压 0.746  1.235

n      = 3;         % 谐波腔次数
n2     = n^2;
k_fp     = sqrt(1/n2-1/(n2-1)*(U0/V_mc)^2);
fais_fp  = pi-asin(n2/(n2-1)*U0/V_mc);
nfaih_fp = atan(tan(fais_fp)/n);

% psih detuning angle
fre_mc = cspeed/C*h;
fre_hc = fre_mc*n;
psih_fp = pi/2 - nfaih_fp;
% loaded shunt impedance
R_L=k_fp*V_mc/(-2*Ib*cos(psih_fp));
disp(['最优负载阻抗 ',num2str(R_L/1e6),' MOhm']);
% assuming the loaded quality factor Q_L = R_L * 80
Q_L=R_L/45;
fre_det=tan(pi-psih_fp)/(2*Q_L)*fre_hc;
disp(['最优失谐频率 ',num2str(fre_det/1e3),' kHz']);

%% near optimum lengthening
Q=5e5;Rs=Q*90;
cos_psih_no=k_fp*V_mc/(-2*Ib*Rs);
delta_fr = -tan(acos(cos_psih_no))*fre_hc/(2*Q);

disp(['近最优失谐频率 ',num2str(delta_fr/1e3),' kHz']);
