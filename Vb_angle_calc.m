function [Vb_angle]=Vb_angle_calc(Re,Im)
% ����г��ǻѹʸ������λ
Vb_angle = atan(Re./Im);
Vb_angle(Im<0) = Vb_angle(Im<0)-pi/2;
Vb_angle(Im>0) = Vb_angle(Im>0)+pi/2;
end