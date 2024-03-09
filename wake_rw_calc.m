function [Wake_rw,wake_bin]=wake_rw_calc(t,bin_tau,tau_q,wakelong)

delta_t = 1e-17;
N = 1e-12/bin_tau;
wake_bin=zeros(N+1,1);
for i = 1:N
timei   = (i-1)*bin_tau+t(1):delta_t:bin_tau*i+t(1);
% 第一个bin 和 第二个bin 用平均值表示
wake_bin(i+1) = sum(interp1(t,wakelong,timei)*delta_t)/bin_tau;
end
wake_bin(1) = 0.5 * wake_bin(2);

Wake_rw = interp1(t,wakelong,tau_q(N+2:end));
Wake_rw =[wake_bin;Wake_rw];
end