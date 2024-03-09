function [Vg_mc_track,fai_s]=Vg_feedback(Vrf_ideal,V_load_cpu_mean,Vg_mc_track,mode)
% test feedback
Vg_ideal = Vrf_ideal-V_load_cpu_mean;
V_condition = abs((Vg_ideal-Vg_mc_track)/Vg_mc_track);
if V_condition > 1e-9
    switch mode
        case 0
            Vg_mc_track = -6.2261e+05 - 1.5251e+05i; % 为设定值，不变值
        case 1
            Vg_mc_track = Vg_ideal;  % is efficient to stablize the beam
        case 2
            Vg_mc_track = Vg_ideal + (Vg_mc_track-Vg_ideal)/10;
        case 3 
            Vg_mc_track = Vg_ideal + (Vg_mc_track-Vg_ideal)/20;
    end 
end
[Vg_angle]       = round(Vb_angle_calc(real(Vg_mc_track),imag(Vg_mc_track))*1e9)/1e9;
fai_s            = pi/2-Vg_angle;
end