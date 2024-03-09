% 多粒子多束团纵向追踪-STABLE
% 研究谐波腔致束团拉伸
% GPU 加速
% 作者： 何天龙
% 时间： 20200706
% time:  20200910 (modified: w/o high Q approximation)
% time : 20210213 (modified: include the BL of main cavity)
% time : 20211121 (modified: include HOMs)
% time : 20220212 (modified: include realistic PI feedback,test)
% time : 20230212 (modified: include realistic PI feedback)
clc;clear;

fre_shift_scan=[60:10:100,150:50:250]*1e3; % scan detuning

for Ii = 1:length(fre_shift_scan)
%% beam parameters
% HALF 参数
cspeed = 299792458; 
sigma_t0 = 10e-12;      % s  束团归一化长度  （用于计算上归一化）
sigma_e0 = 7.44e-4;     %    束团归一化能散 6.598e-4;  1e-3   1.0e-3
alpha_c  = 9.4e-5;      % 动量紧缩因子 6.3e-5; 8.1e-5;9.4e-5
tau_s    = 14e-3;       % 阻尼时间 s 可调  取较小的值可加速收敛 14.2
tau_z    = 14e-3;
I0     = 350e-3;         % 流强 A   350e-3
E0     = 2.2e9;          % 能量 eV
U0     = 400e3;        % 能损 eV  600e3   218e3  198.8 186.5
V_mc   = 1.2e6;          % 主腔腔压 V 1.235e6  0.746e6  0.85  1.44
h      = 800;            % 谐波数
n_hc   = 3;
Q_hc   = 2e8;            % 品质因素 21e3 48.8e3  5e5
R_hc   = Q_hc*39;        % 特征阻抗 200  39,45
C      = 479.86;
fre_shift = fre_shift_scan(Ii);     % 失谐频率 Hz  128300 205940 124300
% fre_shift = detune_HC_calc(I0,n_hc,C,h,U0,V_mc,R_hc,Q_hc);% in near-optimum lengthening condition
% Normal-MC
% Q_mc   = 6095;            % shunt impedance of main cavity
% R_mc   = 6095*119*3;      % quality factor of main cavity
% fre_shift_mc = -54.0e3;   % Hz 主腔失谐频率 根据负载角计算失谐频率
% Super-MC
Q_mc   = 1.1e5;         % shunt impedance of main cavity
R_mc   = Q_mc*45;        % quality factor of main cavity
fre_shift_mc = -7.0e3;   % Hz 主腔失谐频率 根据负载角计算失谐频率  -1.2e3 -2.1

% fill pattern
% pattern(1:10:h)=1;

% hybrid mode
% pattern = zeros(1,h);pattern(1:639)=1;pattern(720)=1;

% standard mode
% pattern  = ones(1,h);
% pat      = 41:50;
% for i = 2:16
%     pat  = [pat,50*i-9:50*i];
% end
% pattern(pat)=0;
% 

% 80% filling rate
pattern  = ones(1,h);
pat      = 33:40;
for i = 2:20
    pat  = [pat,40*i-7:40*i];
end
pattern(pat)=0;

% 90% filling rate
% pattern  = ones(1,h);
% pat      = 37:40;
% for i = 2:20
%     pat  = [pat,40*i-3:40*i];
% end
% pattern(pat)=0;

fillrate = length(find(pattern==1))/h;
HALF = machine(C,I0,U0,E0,tau_s,tau_z,sigma_t0,sigma_e0,alpha_c,h,V_mc,n_hc,R_hc,Q_hc,fillrate,fre_shift,Q_mc,R_mc,fre_shift_mc);
HALF.ShortRange_on = 1; % 0 - neglecting short range effect, 1 considering.

HOMs_m0;  % add HOMs  see the codes of HOM_m0.m
PI_Set;   % add PI    see the codes of PI_Set.m

%% bunch generation
Par_num = 1e4; Bun_num = length(find(pattern==1));

%% charge pattern

% HALF
charge = ones(1,h).*pattern; 

% % charge(720)=5; % hybrid mode, large bunch charge / other bunch charge =5/1
% % charge = charge+TruncatedGaussian(1, [-3,3], [1,h])*0.06.*pattern; % 6% error
charge = charge/sum(charge)*Bun_num;

% generation of initial distribution
if Ii==1
q1 = TruncatedGaussian(1, [-3,3], [Par_num,1]);
p1 = TruncatedGaussian(1, [-3,3], [Par_num,1]);
q  = repmat(q1,1,Bun_num);
p  = repmat(p1,1,Bun_num);%
% CPU to GPU
Q=gpuArray(single(q)); P=gpuArray(single(p));     % single type
end
index_add = 1:Bun_num;
index_add = gpuArray(single(index_add-1));

Dq = 0.05;
%% wake data (intrabunch motion)
tau_q = (0:Dq:200)'*sigma_t0;
Wake_inter_hc = -HALF.wr_hc *  R_hc /Q_hc*exp(-tau_q*HALF.wr_hc/2/Q_hc) .*(...
    cos(tau_q*HALF.wr_hc*HALF.rot_coef_hc)-HALF.VbImagFactor_hc*sin(tau_q*HALF.wr_hc*HALF.rot_coef_hc));
Wake_inter_hc(1) = Wake_inter_hc(1)/2;

Wake_inter_mc = -HALF.wr_mc *  R_mc /Q_mc*exp(-tau_q*HALF.wr_mc/2/Q_mc) .*(...
    cos(tau_q*HALF.wr_mc*HALF.rot_coef_mc)-HALF.VbImagFactor_mc*sin(tau_q*HALF.wr_mc*HALF.rot_coef_mc));
Wake_inter_mc(1) = Wake_inter_mc(1)/2;

Wake_inter = Wake_inter_hc + Wake_inter_mc;

% HOMs
if ~isempty(Q_hom_m0)
for j=1:length(HALF.Q_hom_m0)
    Wake_inter_hom = -HALF.wrf_hom_m0(j)*HALF.R_hom_m0(j)/HALF.Q_hom_m0(j)*...
        exp(-tau_q*HALF.wrf_hom_m0(j)/2/HALF.Q_hom_m0(j)).*(cos(tau_q*...
        HALF.wrf_hom_m0(j)*HALF.rot_coef_hom_m0(j))-HALF.VbImagFactor_hom_m0(j)*...
        sin(tau_q*HALF.wrf_hom_m0(j)*HALF.rot_coef_hom_m0(j)));
    Wake_inter_hom(1) = Wake_inter_hom(1)/2;
    Wake_inter = Wake_inter + Wake_inter_hom;
end
end
bin_tau = Dq*sigma_t0;
%% Longitudinal RW wake
% load('wakez_rw.mat');
% Wake_rw = wake_rw_calc(t,bin_tau,tau_q,wakelong);
% % plot(t,wakelong,'r','LineWidth',1.5);hold on;
% % plot(tau_q,Wake_rw,'b','LineWidth',1.5);xlabel('t [s]');ylabel('V/C');grid minor;title('long.rw wake');
% Wake_inter = Wake_inter - Wake_rw;  % add longitudinal rw wake

%% Longitudinal Geometry wake
% load('wakez_geo.mat');
% Wake_geo = interp1(t,wakelong,tau_q(2:end)); Wake_geo=[wakelong(1)/2;Wake_geo];
% %plot(tau_q,Wake_geo,'b','LineWidth',1.5);xlabel('t [s]');ylabel('wake V/C');grid minor;title('long.geo. wake');
% Wake_inter = Wake_inter - Wake_geo;
%% BBR WAKE
% fr = 30e9; Rs = 2e3;
% Wake_BBR = 2*pi*fr*Rs*exp(-tau_q*2*pi*fr/2).*(cos(tau_q*2*pi*fr*sqrt(0.75))-sin(tau_q*2*pi*fr*sqrt(0.75))/2/sqrt(0.75));
% Wake_BBR(1) = Wake_BBR(1)/2;
% % plot(tau_q,Wake_BBR);
% Wake_inter = Wake_inter - Wake_BBR;

% wr_bb = 11.549e9*2*pi;   % Hz
% Rs_bb = 5.730e3;         % Ohm
% Q_bb  =6;
% az = wr_bb/(2*Q_bb);wr1 = sqrt(wr_bb^2-az^2);
% Wake_bb = -2*az*Rs_bb*exp(-az*tau_q).*(cos(-wr1*tau_q)+az/wr1*sin(-wr1*tau_q));
% Wake_bb(1) = Wake_bb(1)/2;
% Wake_inter = Wake_inter + Wake_bb;
%%
Wake_inter    = gpuArray(single(Wake_inter));
% plot(tau_q,Wake_inter);

%% start tracking Track_num = 1e3
% charge per macro-particle   : HALF.qc
% 不等电荷量填充时，每个束团的宏粒子电荷量不等，注意区别
HALF.qc   = charge.*pattern * HALF.qc / Par_num;              %由单个元素变为一行矩阵
% induced voltage per macro-particle  : HALF.V_b
HALF.Vb_hc  = HALF.qc * HALF.wr_hc * HALF.R_hc / HALF.Q_hc; 
HALF.Vb_mc  = HALF.qc * HALF.wr_mc * HALF.R_mc / HALF.Q_mc; 

% HALF.V_b  = HALF.qc * HALF.w_r * HALF.R_hc / HALF.Q_hc *(1+1i*HALF.VbImagFactor); 
% initial loaded voltage
if Ii==1
V_hc_load_0_real = real(HALF.V_hc_load_0);
V_hc_load_0_imag = imag(HALF.V_hc_load_0);
V_hc_load_0      = V_hc_load_0_real+1i*V_hc_load_0_imag;
% V_hc_load_0      =0;

V_mc_load_0_real = real(HALF.V_mc_load_0);
V_mc_load_0_imag = imag(HALF.V_mc_load_0);
V_mc_load_0      = V_mc_load_0_real+1i*V_mc_load_0_imag;
% V_mc_load_0=0; 
end

rot_decay_coef_hc = 1i * HALF.rot_coef_hc - 1 / (2 * HALF.Q_hc);  % 旋转 衰减项
TbAng_coef_hc     = exp(rot_decay_coef_hc * HALF.angle_hc);       % 注意符号正负
exp_ang_coef_hc   = -rot_decay_coef_hc * HALF.wr_hc * sigma_t0;

rot_decay_coef_mc = 1i * HALF.rot_coef_mc - 1 / (2 * HALF.Q_mc);  % 旋转 衰减项
TbAng_coef_mc     = exp(rot_decay_coef_mc * HALF.angle_mc);       % 注意符号正负
exp_ang_coef_mc   = -rot_decay_coef_mc * HALF.wr_mc * sigma_t0;

if ~isempty(Q_hom_m0)
HALF.Vb_hom  = HALF.qc .* (HALF.wrf_hom_m0 .* HALF.R_hom_m0 ./ HALF.Q_hom_m0); 
rot_decay_coef_hom = 1i * HALF.rot_coef_hom_m0 - 1 ./ (2 * HALF.Q_hom_m0);   % 旋转 衰减项
TbAng_coef_hom     = exp(rot_decay_coef_hom .* HALF.angle_hom_m0);           % 注意符号正负
exp_ang_coef_hom   = -rot_decay_coef_hom .* HALF.wrf_hom_m0 * sigma_t0;
V_hom_load_0_real = zeros(HALF.Q_hom_length,1);
V_hom_load_0_imag = zeros(HALF.Q_hom_length,1);
V_hom_load_0      = V_hom_load_0_real+1i*V_hom_load_0_imag;
end

wake_kick_coef = HALF.qc * HALF.kick_coef;

% 发射机电压矢量替代之前的Vrf矢量
Vg_mc = abs(HALF.Vg_mc_init);
[Vg_angle]=round(Vb_angle_calc(real(HALF.Vg_mc_init),imag(HALF.Vg_mc_init))*1e12)/1e12;
HALF.Vg_mc_track = HALF.Vg_mc_init;
HALF.rfcoef1_track     = HALF.rfcoef1 / HALF.V_mc * Vg_mc;
fai_s_track            = pi/2-Vg_angle;                    % 发射机电压矢量的同步相位

Track_num  = 10e4;   % set tracking turns
% record parameters 
Recor_step = 10;
HALF.Recor_step=Recor_step;
Recor_num  = Track_num / Recor_step;
% Vg_mc_track_record = zeros(1,Track_num*h/10);
Vb_hc_track_record = zeros(1,Track_num);
record_Q_mean = zeros(Recor_num,Bun_num);record_Q_std = zeros(Recor_num,Bun_num);
record_P_mean = zeros(Recor_num,Bun_num);record_P_std = zeros(Recor_num,Bun_num);

%%
gd = gpuDevice(); 
tic;

pikp = [1];   % kp  of PI control
piki = [1e-5];% ki  of PI control
j=1;
for pikpi = 1:length(pikp)
PI.KP = pikp(pikpi);
for pikii = 1:length(piki)
record_th = 0;
PI.KI = piki(pikii);
HALF.Ig_track=zeros(1,HALF.h/PI.m*Track_num);
Ig_track_tnum = 0;
for i =1:Track_num
    % drift
    Q = Q + P * HALF.drift_coef;   
    Q_min = min(Q);    
    Q_new = round((Q - Q_min) *(1/Dq));
        
%%  HOMs Kick
    if ~isempty(Q_hom_m0)
        for jj=1:HALF.Q_hom_length
            exp_angle=exp(exp_ang_coef_hom(jj) * Q);
            exp_angle_sum=gather(sum(exp_angle));
            V_load_cpu = double(exp_angle_sum.*HALF.Vb_hom(jj,pattern==1));% *HALF.Vb_hom
            [V_load,V_hom_load_0(jj)]=VoltageLoadCalc_matlab(V_hom_load_0(jj),V_load_cpu,TbAng_coef_hom(jj),pattern);                

            V_load_cpu = V_load(pattern==1) * HALF.kick_coef;   % 约化V_load;
            V_load     = gpuArray(single(V_load_cpu)); 
            V_hom_load_kick = V_load./exp_angle;
            P = P - real(V_hom_load_kick) + imag(V_hom_load_kick) * HALF.VbImagFactor_hom_m0(jj);
        end
    end
%% Harmonic cavity     
    % beam induced voltage at nominal bucket position HHC
    exp_angle  = exp(exp_ang_coef_hc * Q);
    exp_angle_sum= gather(sum(exp_angle));       % 耗时 0.007s sum()函数较慢    
    V_load_cpu = double(exp_angle_sum.*HALF.Vb_hc(pattern==1)); % *HALF.Vb_hc  
    [V_load,V_hc_load_0]=VoltageLoadCalc_matlab(V_hc_load_0,V_load_cpu,TbAng_coef_hc,pattern); 
    Vb_hc_track_record(i)=mean(V_load);
%  谐波腔腔压矢量图示
    if mod(i,1000)==0
        figure(13)
        subplot(2,1,1)
        plot(abs(V_load)/1e3);
        title('Harmonic Cavity');ylabel('Amplitude [kV]');
        subplot(2,1,2)
        plot(angle(V_load)/pi*180);ylabel('Phase [deg]');xlabel('Bucket ID');
    end
    V_load_cpu = V_load(pattern==1)*HALF.kick_coef;   % 约化V_load;
    V_load     = gpuArray(single(V_load_cpu));    
    % intrabunch kick    - V_load_kick    real part
    V_hc_load_kick = V_load./exp_angle;   
% _________________________________________________________________________    
%% Main cavity       
    % beam induced voltage at nominal bucket position MC
    exp_angle  = exp(exp_ang_coef_mc * Q);
    exp_angle_sum= gather(sum(exp_angle));       % 耗时 0.007s sum()函数较慢    
    V_load_cpu = double(exp_angle_sum.*HALF.Vb_mc(pattern==1)); % *HALF.Vb_mc
    [Vc_mc,Vg_mc_track,HALF.Vg_mc_track_0,V_load,HALF.V_mc_load_0,PI]=PI_Control(PI,HALF.Vrf_ideal,HALF.Vg_mc_track_0,...
    HALF.V_mc_load_0,V_load_cpu,TbAng_coef_mc,pattern); % every 5120 buckets to do PI

%     if i == 3e4
%         HALF.Vrf_ideal=HALF.Vrf_ideal*1.08;   % 测试PI反馈对腔压设定值响应能力
%     end

    Ig_track_num = length(PI.Ig_track);
    if Ig_track_num>10000
        Ig_track_tnum = Ig_track_tnum + Ig_track_num;
        HALF.Ig_track(Ig_track_tnum-Ig_track_num+1:Ig_track_tnum)=PI.Ig_track;
        PI.Ig_track = [];
    end

    V_mc_kick = gpuArray(single(Vc_mc(pattern==1)*HALF.kick_coef))./exp_angle;
%   主腔腔压矢量图示    
    if mod(i,1000)==0
        figure(15)
        subplot(2,1,1)
        plot(abs(Vc_mc)/1e3);title('Main Cavity');ylabel('Amplitude [kV]');
        subplot(2,1,2)
        plot(angle(Vc_mc)/pi*180);ylabel('Phase [deg]');xlabel('Bucket ID');

        % generator current
        figure(666);
        subplot(2,1,1);plot(abs(PI.Ig_track));ylabel('amplitude');title('generator current');
        subplot(2,1,2);plot(angle(PI.Ig_track));ylabel('phase');
        % 发射机功率
        % Pg_mc = 1/8*PI.Ig_track.^2*R_mc_0/betacoupling*4;
%         figure(667);plot(abs(Pg_mc)/1e3);ylabel('P_g  [kW]');
    end
%% short-range wake kick   
    % count bins
    if HALF.ShortRange_on ==1                 % modified in 2022/11/14
        binnum=max(max(Q_new))+1; binnum=gather(binnum);
        bin_num_q=sum(BinNumCalZ(binnum,Q_new));   % double type
        bin_num_q=reshape(bin_num_q,binnum,Bun_num);    
        kick_conv = conv2(bin_num_q,Wake_inter(1:binnum));    
        kick_conv = kick_conv(1:binnum,:) .* wake_kick_coef(pattern==1)*min(i/5000,1);
    % wake kick
        Q_new = Q_new + (1 + index_add * binnum);  % modified in 2020/09/22
        wake_kick = kick_conv(Q_new);
    else
        wake_kick = 0;
    end
    % radiation damping and quantum excitation term  + wake_kick
    rad_quan_kick = -HALF.radampcoef * P + HALF.quanexcoef *...
        gpuArray.randn(Par_num,Bun_num,'single');

    P = P + rad_quan_kick - HALF.ploss;

    P = P - real(V_hc_load_kick) + imag(V_hc_load_kick) * HALF.VbImagFactor_hc...
        - real(V_mc_kick) + imag(V_mc_kick) * HALF.VbImagFactor_mc+ wake_kick;    
    
    if mod(i,2000)==0
        Centroid_std=std(record_Q_mean(record_th,:))*HALF.sigma_t0*1e12;
        disp(['tracking turn = ',num2str(i),'; Centroid_std = ',num2str(Centroid_std),'ps']);
        toc;
    end
    % output data
    if mod(i,Recor_step)==0
        record_th = record_th +1;
        record_Q_mean(record_th,:)=gather(mean(Q));
        record_Q_std(record_th,:)=gather(std(Q));
        record_P_mean(record_th,:)=gather(mean(P));
        record_P_std(record_th,:)=gather(std(P));
    end
end
wait(gd);
toc;
%%
filename=['HALF_80percent_I0',num2str(I0*1e3),'mA','_RLfp',num2str(R_hc),...
    '_QLfp',num2str(Q_hc),'_detune',num2str(fre_shift),'_kp',num2str(PI.KP),'_ki',num2str(PI.KI),'_',num2str(1),'.mat'];
save(filename,'record_Q_mean','record_Q_std','record_P_mean','record_P_std','Q','Track_num','HALF','PI','Bun_num','Vb_hc_track_record');
end
end
end
%%
figure(1);
Recor_step=HALF.Recor_step;
Nturns = (1:Track_num/Recor_step)*Recor_step;
for i=1:2:8
    subplot(2,2,1)
    plot(Nturns,record_Q_mean(:,i)*HALF.sigma_t0*1e12); hold on;
    subplot(2,2,2)
    plot(Nturns,record_Q_std(:,i)*HALF.sigma_t0*1e12); hold on;
    subplot(2,2,3)
    plot(Nturns,record_P_mean(:,i)*HALF.sigma_e0); hold on;
    subplot(2,2,4)
    plot(Nturns,record_P_std(:,i)*HALF.sigma_e0); hold on; 
end
subplot(2,2,1);ylabel('<\tau>  [ps]');xlabel('turns');xlim([1,Track_num]);grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
subplot(2,2,2);ylabel('\sigma_{\tau}  [ps]');xlabel('turns');xlim([1,Track_num]);grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
subplot(2,2,3);ylabel('<\delta> ');xlabel('turns');xlim([1,Track_num]);grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
subplot(2,2,4);ylabel('\sigma_{\delta} ');xlabel('turns');xlim([1,Track_num]);grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
%%
% 统计沿着束团 长度分布，中心分布  1:100:2000
figure(2);
% for i=10
% subplot(2,1,2);plot(mean(record_Q_mean(end-i:end,:))*HALF.sigma_t0*1e12,'.');hold on;
% ylabel('<\tau>  [ps]');xlabel('bunch number');
% subplot(2,1,1);plot(mean(record_Q_std(end-i:end,:))*HALF.sigma_t0*1e12,'.');hold on;
% ylabel('\sigma_{\tau}  [ps]');xlabel('bunch number');
% end
for i=0
subplot(1,2,2);plot(record_Q_mean(end-i,:)*HALF.sigma_t0*1e12,'.');hold on;
ylabel('<\tau>  [ps]');xlabel('bunch number');
subplot(1,2,1);plot(record_Q_std(end-i,:)*HALF.sigma_t0*1e12,'.');hold on;
ylabel('\sigma_{\tau}  [ps]');xlabel('bunch number');
end
mean(record_Q_std(end-i,:)*HALF.sigma_t0*1e12)
subplot(1,2,2);
% ylim([-15,15]);
% grid minor;
set(gca,'FontName','Times New Roman','FontSize',12);xlim([1,Bun_num]);
subplot(1,2,1);
% grid minor;
set(gca,'FontName','Times New Roman','FontSize',12);xlim([1,Bun_num]);
%%  统计作密度分布图  Dq = 0.4;
Dq = 0.5;
Q_min = min(Q); 
tau_min = gather(Q_min)*HALF.sigma_t0;
Q_new = round((Q - Q_min) *(1/Dq));
binnum=max(max(Q_new))+1; binnum=gather(binnum);
bin_num_q=sum(BinNumCalZ(binnum,Q_new));
bin_num_q=reshape(bin_num_q,binnum,Bun_num);
figure(4);
bin_i = [1];
colorset =[1 0 0;0 1 0;0 0 1;0 0 0.5;1 0.5 0.5; 0.5 1 0.5;0.5 0.5 1;1 0 1;0 1 1;1 1 0];
for i =1:length(bin_i)
    bin_range = (1:binnum)*Dq*HALF.sigma_t0+tau_min(bin_i(i));
    plot((bin_range'-Dq*HALF.sigma_t0)*1e12,bin_num_q(:,bin_i(i))/1e4/5,'-','Color',colorset(i,:),'Linewidth',2);hold on; % rs
end
ylabel('norm.density ');xlabel('\tau [ps]'); 
xlim([-150,150]);
set(gca,'FontName','Times New Roman','FontSize',14);
%% 画出Vg电压
figure(656)
plot(abs(Vb_hc_track_record)/1e3/mean(abs(Vb_hc_track_record(1:50000))/1e3));title('Harmonic Cavity');hold on;
plot(angle(Vb_hc_track_record)/pi*180/mean(angle(Vb_hc_track_record(1:50000))/pi*180));hold on;
plot(Nturns,record_P_mean(:,1)/7+1);hold on;
legend('Amplitude','Phase','<\delta>');
xlabel('Turns');ylabel('Norm.Amp. [a.u.]');xlim([0,10e4]);
%%
figure(666)
plot(Vb_hc_track_record(10000:50000)-1i*mean(imag(Vb_hc_track_record(10000:100000))));hold on;
xlabel('Real part [V]');
ylabel('Imag part [V]');
plot(Vb_hc_track_record(50000:100000)-1i*mean(imag(Vb_hc_track_record(50000:100000))));hold on;
%% 主腔发射机功率  PI.Ig_track 
% 发射机电流
figure(666)
subplot(2,1,1);
plot(abs(PI.Ig_track));hold on
subplot(2,1,2);
plot(angle(PI.Ig_track));hold on;
%% 发射机功率
Q_mc_0 = 5e8;R_mc_0 = Q_mc_0*44.5; betacoupling = Q_mc_0/HALF.Q_mc-1;% main cavity param.
figure(667)
Pg_mc = 1/8*HALF.Ig_track.^2*R_mc_0/betacoupling*4; % *4 due to similar to Ib
plot(abs(Pg_mc)/1e3);hold on;ylabel('P_g  [kW]');

%% FFT分析振荡频率 Q
% for i=1:1
% mean_q = record_Q_mean(5000*(i)+1:5000*(i+1),1)'; % 1 first bunch
% mean_q = mean_q-mean(mean_q); % 去DC
% n_turns= length(mean_q);
% % 统计质心的振荡频率
% % 注意此处是每10圈记录一次数据
% freqs = 0.00001:1/n_turns:0.5;amp = abs(fft(mean_q));
% figure(12)
% % plot(freqs/10,amp(1:length(freqs))/max(amp(1:length(freqs)))); %/10 表示每10圈记录一次数据
% % xlim([0,0.5]);
% plot((freqs/2)*(299792458/HALF.C),amp(1:length(freqs))/max(amp(1:length(freqs)))); %/10 表示每10圈记录一次数据
% % xlim([25e3,35e3]);
% hold on;
% pause(1)
% end
