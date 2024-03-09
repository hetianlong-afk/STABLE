% fft analysis
% FFT分析频散
% 分析纵向振荡频率及其分散
% 分析Ig调制频率
% 需要统计每圈数据，圈数是2的指数
%% centroid
clc;clear;
kp_scan = [2.2];
for kpi = 1:length(kp_scan)
kp =  kp_scan(kpi);
filename=['.\HALF_Hybrid_I0350mA_RLfp19500000_QLfp500000_detune40000_kp',num2str(kp),'_ki1e-06.mat'];
load(filename);
mean_q = record_P_mean(2000:5000,1)'; % 1 first bunch
n_turns= length(mean_q);
% 统计质心的振荡频率
% 注意此处是每10圈记录一次数据
freqs = 0.00001:1/n_turns:0.5;amp = abs(fft(mean_q));
figure(1)
plot(1./(freqs/10),amp(1:length(freqs))); %/10 表示每10圈记录一次数据
xlim([0,2000]);
hold on;
end
%% Ig
Q_mc_0 = 2e9;R_mc_0 = Q_mc_0*45; betacoupling = Q_mc_0/HALF.Q_mc-1;% main cavity param.
Pg_mc = 1/8*HALF.Ig_track.^2*R_mc_0/betacoupling*4; % *4 due to similar to Ib
Pg_mcabs = abs(Pg_mc)/1e3;
Pg_mcabs_mean = Pg_mcabs-mean(Pg_mcabs);
figure(2)
plot(abs(Pg_mc)/1e3);hold on;ylabel('P_g  [kW]');
