function loading_angle_calc(Vg_mc_track,V_load_cpu_mean,det_angle_mc)
% 计算负载角
% Vg_mc_track ： 发射机电压矢量
% V_load_cpu_mean： 各束团平均负载电压矢量
V_mc_total = Vg_mc_track+V_load_cpu_mean;

angle_V_mc = Vb_angle_calc(real(V_mc_total),imag(V_mc_total));
angle_V_ge = Vb_angle_calc(real(Vg_mc_track),imag(Vg_mc_track));
fais_mc_whc = pi/2 - angle_V_mc;

angle_load = fais_mc_whc-(pi/2 - angle_V_ge-det_angle_mc);
disp(['初始负载角设置为：',num2str(angle_load/pi*180),'deg']);
end