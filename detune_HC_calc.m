function [fre_detRQ] = detune_HC_calc(Ib,n_hc,C,h,U0,V_mc,R_hc,Q_hc)
n2          = n_hc^2;
fre_hc      = 299792458/C*h*n_hc; % 谐波腔的谐振频率 Hz
k_fp        = sqrt(1/n2-1/(n2-1)*(U0/V_mc)^2);
phih_RQ     = acos(k_fp*V_mc./(-2*Ib*R_hc));
fre_detRQ   = round(tan(pi-phih_RQ)./(2*Q_hc)*fre_hc);
disp([num2str(Ib*1e3), 'mA', '对应的失谐频率：',num2str(fre_detRQ),'Hz']);
end
