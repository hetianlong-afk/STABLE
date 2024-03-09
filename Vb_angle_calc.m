function [Vb_angle]=Vb_angle_calc(Re,Im)
% 计算谐波腔压矢量的相位
Vb_angle = atan(Re./Im);
Vb_angle(Im<0) = Vb_angle(Im<0)-pi/2;
Vb_angle(Im>0) = Vb_angle(Im>0)+pi/2;
end