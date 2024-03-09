%% 计算增长率
turn = 50000:100000;
x=record_P_mean(turn,1)';
pass = turn*2;

xx=log(abs(x));
TF = islocalmax(xx);

% figure(50)
% plot(pass,xx,pass(TF),xx(TF));

a=polyfit(pass(TF),xx(TF),1);

xs = exp(a(2)+a(1)*pass);

figure(100)
plot(pass,x,pass,xs,'.');xlabel('Turns');ylabel('<\delta>');
legend('Tracking','Fitted');
disp(['The growth rate is ',num2str(a(1)/HALF.T0)]);
disp(['The radiation damping rate is ',num2str(1/HALF.tau_z)]);
