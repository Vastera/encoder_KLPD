[phase,speed]=simulated_encoder_sig;
fs=2e4;t=1/fs:1/fs:1;
figure('Name','Generated Speed')
plot(t,speed);ylim([270,330]);
figure('Name','Generated phase')
plot(t,phase);
%%%% test KLPD %%%%%%%%%%%
M=1:400;
[M_opt,kurt]=KLPD(phase,fs,max(M));
figure('Name','Kurtosis for different filter length')
plot(M*2+1,kurt);
hold on;plot(M_opt*2+1,kurt(M_opt),'ro');
text(M_opt*2+1,kurt(M_opt)+1,[num2str(M_opt*2+1),' , ',num2str(kurt(M_opt)+1)]);
xlabel('Filter length');ylabel('Kurtosis');

%%%%%%%%%%% use the optimal length %%%%%%%%%%%%%
figure('Name','fltered speed')
v=LPD(phase,M_opt,fs);
plot(t(1:length(v)),circshift(v,50));hold on; plot(t,speed);
xlabel('Time [s]');ylabel('IAS [rpm]');
ylim([270 330]);
SetFigureProperties;