%%% TESTES ESTIMADOR DFT-PRONY

clc
clear all
close all


set(0, 'defaultAxesFontSize',12);
set(0, 'defaultAxesFontName','times');




% ---------------------------------------------------------------------
% INICIALIZAÇÕES

N = 256;       % número de amostras
m = 128;       % passo de janelamento
Fs = 6400;    % Taxa de amostragem
f0 = 50;      % freq. fundamental (teórica)
noise = 60;    % nível de ruído em dB ( 0 = sem ruído)
T = 1.2;      % tempo total (em segundos)
param = 0;    % amort(w);    % freq. de modulação, por exemplo
type = 0;     % seleciona tipo de teste (ver gerador_sinais.m)

%f1 = 50.12;
%freqs = [50.23];        % Vetor de frequências
freqs = [50. 150. 251. ];        % Vetor de frequências
%freqs = [f1 3*f1 5*f1 7*f1 9*f1]

phases(1,:) = [3.96 4.85 4.75 3.11 2.1];
phases(2,:) = phases(1,:)-2*pi/3;
phases(3,:) = phases(1,:)+2*pi/3;
%phases = zeros(3,2);
% phases = [0 0 0 0;           % vetor de fases
%          -2*pi/3 0 0 0;
%          +2*pi/3 0 0 0];
% phases(1,:) = [3.96 4.85 4.75 3.11 6.01 5.54 0.16 3.91];           % vetor de fases
% phases(2,:) = phases(1,:)-2*pi/3;
% phases(3,:) = phases(1,:)+2*pi/3;
% phases = [5.1 0 0 0 6.01 5.54 4.16 3.91;  
%          -2*pi/3 0 0 0 0 0 0 0;
%          +2*pi/3 0 0 0 0 0 0 0];    
% phases(1,1) = deg2rad(20);
% phases(2,1) = deg2rad(-100);
% phases(3,1) = deg2rad(+140);
% 
% phases(1,2) = deg2rad(0);
% phases(2,2) = deg2rad(-120);
% phases(3,2) = deg2rad(+120);

% amps(1,:) = [1 0.5 0.3 0.2 0.1];
% amps(2,:) = amps(1,:);
% amps(3,:) = amps(1,:);
%amps = [1 0.1; 1 0.1; 1 0.1;];
% amps = [6 2.21 0.812 0.299];           % vetor de amplitudes
% amps(2,:) = amps(1,:);
% amps(3,:) = amps(1,:);

amps(1,:) = [1 0.3 0.1];           % vetor de amplitudes
amps(2,:) = amps(1,:);
amps(3,:) = amps(1,:);
%         1 0.3 0 0 0 0 0 0;
%         1 0.3 0 0 0 0 0 0]; 

% ------------------------------------------------------------------------------------------------------------------
% GERAÇÃO DO SINAL

[Va, Vb, Vc, t, refs, mod] = gerador_sinais(f0, freqs, phases, amps,Fs, T,noise, type, param);

% figure, hold on, grid on
% % plot(t,mod*10)
% plot(t,Va)
% plot(t,Vb)
% plot(t,Vc)
% legend('Va','Vb','Vc'), title('Sinais')
% xlabel('tempo (s)'), ylabel('amplitude')

%-------------------------------------------------------------------------------------------------------------------
% SEGMENTAÇÃO (JANELAMENTO)
[Va_seg, Vb_seg, Vc_seg, t_seg, ref_seg] = segmentador(Va, Vb, Vc, t, refs, N, m);



% Plot amostragem (ciclo a ciclo)
w = 4;   % numero de janelas no plot
figure
for i=1:w
    subplot(w,1,i), hold on, grid on
    plot(t_seg(i,:),Va_seg(i,:),'o-')
end

% figure %Plot referências
% subplot(3,1,1), hold on, grid on
% plot(t_seg(:,1), ref_seg(1,:), 'o-')
% xlabel('tempo (s)'), ylabel('freq (Hz)'), title('Referências (segmentada)')
% subplot(3,1,2), hold on, grid on
% plot(t_seg(:,1),ref_seg(2,:), 'o-')
% plot(t_seg(:,1),ref_seg(3,:), 'o-')
% plot(t_seg(:,1),ref_seg(4,:), 'o-')
% xlabel('tempo (s)'), ylabel('amp (pu)')
% subplot(3,1,3), hold on, grid on
% plot(t_seg(:,1),ref_seg(5,:), 'o-')
% plot(t_seg(:,1),ref_seg(6,:), 'o-')
% plot(t_seg(:,1),ref_seg(7,:), 'o-')
% xlabel('tempo (s)'), ylabel('fase (rad)')

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR FASORIAL DFT
[Sa, Sb, Sc, Aa, Ab, Ac, phia, phib, phic] = estimador_dft(Va_seg, Vb_seg, Vc_seg, Fs, f0);

granN = (Fs/N);
gradeN = 0:granN:granN*(N-1);

% % Plot espectro
% figure
% hold on, grid on
% stem(gradeN,abs(Sa(1,:)))
% % stem(abs(Sa))
% % stem(abs(Sa))
% title('Amplitude (pu) - DFT')
% ylabel('pu'), xlabel('Frequência (Hz)')


% CÁLCULO DE COMPONENTES SIMÉTRICAS

bin = round(f0/granN)+1;       % local da componente fundamental no espectro

% COM DADOS FFT
[A0, A1, A2, phi0, phi1, phi2] = comp_simetricas(Aa(:,bin)', Ab(:,bin)', Ac(:,bin)', phia(:,bin)', phib(:,bin)', phic(:,bin)');

% COM DADOS DE REFERÊNCIA (para teste do estimador)
[A0r, A1r, A2r, phi0r, phi1r, phi2r] = comp_simetricas(ref_seg(2,:), ref_seg(3,:), ref_seg(4,:),ref_seg(5,:), ref_seg(6,:), ref_seg(7,:));


% [S0, S1, S2] = comp_simet_all_freq(Sa, Sb, Sc);
% 
% A0 = 2*abs(S0(:,3))/N;
% A1 = 2*abs(S1(:,3))/N;
% A2 = 2*abs(S2(:,3))/N;
% 
% phi0 = angle(S0(:,3));
% phi1 = angle(S1(:,3));
% phi2 = angle(S2(:,3));


% Plot fasores da componente fundamental
figure
subplot(2,3,1), hold on, grid on, plot(A0r, 's-'), plot(A0, 'o-'), title('Zero (pu)'),ylabel('amp. (pu)'), xlabel('Ciclo')
subplot(2,3,4), hold on, grid on, plot(rad2deg(phi0r), 's-'), plot(rad2deg(phi0), 'o-'), ylabel('fase (º)'), xlabel('Ciclo')
subplot(2,3,2), hold on, grid on, plot(A1r, 's-'), plot(A1, 'o-'), title('Positiva (pu)'), ylabel('amp. (pu)'), xlabel('Ciclo')
subplot(2,3,5), hold on, grid on, plot(rad2deg(phi1r), 's-'), plot(rad2deg(phi1), 'o-'), ylabel('fase (º)'), xlabel('Ciclo')
subplot(2,3,3), hold on, grid on, plot(A2r, 's-'), plot(A2, 'o-'), title('Negativa (pu)'), ylabel('amp. (pu)'), xlabel('Ciclo')
subplot(2,3,6), hold on, grid on, plot(rad2deg(phi2r), 's-'), plot(rad2deg(phi2), 'o-'),  ylabel('fase (º)'), xlabel('Ciclo')


%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR DE FREQUÊNCIA - Derivada do ângulo do fasor de seq. positiva

% COM DADOS DE REFERÊNCIA (teste do estimador)
[freq_final_ref] = estimador_freq_deriv_ang(phi0r, phi1r, phi2r, f0, Fs, ref_seg, m, N);

% COM DADOS FFT
[freq_final_fft] = estimador_freq_deriv_ang(phi0, phi1, phi2, f0, Fs, ref_seg, m, N);

figure, hold on, grid on
plot(ref_seg(1,1:end), 'o-')
plot(freq_final_ref, 'o-')
plot(freq_final_fft, 'o-')
% plot(50.23 - d_phi2*0.01 -0.0005, 'o-')
% plot(50.23 + A2-0.0023, 'o-')
title('Frequência (Hz)'), legend('Referência', 'Estimativa (ref)', 'Estimativa (DFT)');
ylabel('Hz'), xlabel('Ciclo')


%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR ESPARSO - DIC. FREQUÊNCIA - OMP_2

P = 100;
Va_seg1 = Va_seg(1,:);
zdft = fft(Va_seg1,N*P);
amp_zdft = 2*abs(zdft)/(N);
ang_zdft = angle(zdft);


% Ajusta amplitude do espectro fourier
Sa = 2*Sa/N;
K = 6; 

P = 100;
H = N*P;
granH = (Fs/(H)) % Nova granularidade
gradeH = 0:granH:granH*((H)-1);


load('C:\dicionarios_cs_dft\block_dic256u.mat')
[Dic] = dic_freq(N, H);

P = 100;
omp = 4;
[Saf1, Aaf1, phiaf1, indx_af1] = estimador_freq_omp(Sa, K, P, omp, Fs, freq_final_ref, Dic);
%Aaf1 = Aaf1*2;
clear Dic % libera memoria


P = 10;
omp = 4;
[Saf2, Aaf2, phiaf2, indx_af2] = estimador_freq_omp_A(Sa, K, P, omp, Fs, f0, freq_final_ref);

%Aaf2 = Aaf2*2;
load('C:\dicionarios_cs_dft\block_dic256us.mat')
P = 10;
omp = 5;
[Saf3, Aaf3, phiaf3, indx_af3] = estimador_freq_omp_B(Sa, K, P, omp, Fs, f0, freq_final_ref, block);

clear block % libera memoria

%-------------------------------------------------------------------------
% ZERO PADDED DFT
P = 10;
granZ = (Fs/(P*N)) % Nova granularidade
gradeZ = 0:granZ:granZ*((P*N)-1);

Va_seg1 = Va_seg(1,:);
% for i=1:P-1
%     Va_seg1 =[Va_seg1 Va_seg];
% end

zdft = fft(Va_seg1,N*P);
amp_zdft = 2*abs(zdft)/(N);
ang_zdft = angle(zdft);


amp_Sa1 = abs(Saf1(1,:));
ang_Sa1 = angle(Saf1(1,:));
amp_Sa2 = abs(Saf2(1,:));
ang_Sa2 = angle(Saf2(1,:));
amp_Sa3 = abs(Saf3(1,:));
ang_Sa3 = angle(Saf3(1,:));

%a_Sa4 = abs(Saf4);
x = 250;
x1 = 2500;
%x=344*4;
figure, hold on, grid on
plot(gradeZ(1:x),amp_zdft(1:x))
plot(gradeH(1:x1),amp_Sa1(1:x1))
plot(gradeZ(1:x),amp_Sa2(1:x))
plot(gradeZ(1:x),amp_Sa3(1:x))
%plot(gradeH(1:x),a_Sa4(1:x))
ylabel('amplitude (p.u.)'),xlabel('frequência (Hz)')
%legend('Zero-padded DFT','Compressive Sensing')

% %--------------------------------------------------------------------------------------------------------------------
% % CÁLCULO DE ERROS
% 
% % índice da dft, alterar conforme incializações
% hf(1) = 1 + 1; % fundamental
% 
% i=1;  % inicio
% [f,~]=size(Aa); % fim (numero de ciclos)
% 
% % for h=0:0;
% % % amplitudes
% % amp_med(1,:) = Aa(i:f,hf(h+1));    % dft
% % amp_med(2,:) = Aaf1(i:f,h*2+1);
% % amp_med(3,:) = Aaf2(i:f,h*2+1); 
% % amp_med(4,:) = Aaf3(i:f,h*2+1);
% % 
% % % fases
% % fase_med(1,:) = phia(i:f,hf(h+1));   
% % fase_med(2,:) = phiaf1(i:f,h*2+1);
% % fase_med(3,:) = phiaf2(i:f,h*2+1);
% % fase_med(4,:) = phiaf3(i:f,h*2+1);
% % 
% % [tve, amp_error, phase_error, freq_error] = calcula_erros(amp_med, fase_med, freq_final_fft, ref_seg(:,i:f));
% % 
% % end
% % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%% PLOT RESULTADOS
% % % 
% % figure, hold on, grid on
% % plot(t_seg(i:f,1), tve(1,:), '*-')
% % plot(t_seg(i:f,1), tve(2,:), 'x-')
% % plot(t_seg(i:f,1), tve(3,:), 'o-')
% % plot(t_seg(i:f,1), tve(4,:), 's-')
% % legend('DFT', 'CS-DE', 'CS-DA', 'CS-DM')
% % ylabel('TVE (%)'), xlabel('tempo (s)')
% 
% 
% 
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % for i=1:4;
% %    max_tve(i) = max(tve(i,:));
% %    max_amp_error(i) = max(amp_error(i,:)); 
% %    max_phi_error(i) = max(phase_error(i,:)); 
% %    
% %    std_tve(i) = std(tve(i,:));
% %    std_amp_error(i) = std(amp_error(i,:)); 
% %    std_phi_error(i) = std(phase_error(i,:)); 
% % end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%% ANÁLISE DAS HARMÔNICAS
% 
% % índice da dft, alterar conforme incializações
% n_cic = round(f0*N/Fs);
% 
% hf(1) = 1 + n_cic*1; % fundamental
% hf(2) = 1 + n_cic*3; % 3rd harmonica
% hf(3) = 1 + n_cic*5; % 5th harmonica
% hf(4) = 1 + n_cic*7; % 7th harmonica
% hf(5) = 1 + n_cic*9; % 9th harmonica
% 
% % hf(1) = 35;
% % hf(2) = 51;
% % hf(3) = 60;
% % hf(4) = 100;
% % hf(5) = 123;
% % hf(6) = 152;
% % hf(7) = 169;
% % hf(8) = 202;
% % hf(1) = 4;
% % hf(2) = 6;
% % hf(3) = 7;
% % hf(4) = 11;
% % hf(5) = 13;
% % hf(6) = 16;
% % hf(7) = 18;
% % hf(8) = 21;
% 
% i=1;  % inicio
% [f,~]=size(Aa); % fim (numero de ciclos)
% 
% 
% for h=0:K-1;
% 
% % amplitudes
% amp_med(1,:) = Aa(i:f,hf(h+1));    % dft
% %amp_med(2,:) = Aaf1(i:f,h*2+1);
% amp_med(2,:) = (Aaf1(i:f,h*2+1)+Aaf1(i:f,h*2+2))/2;
% amp_med(3,:) = Aaf2(i:f,h*2+1); 
% amp_med(4,:) = Aaf3(i:f,h*2+1);
% 
% % fases
% fase_med(1,:) = phia(i:f,hf(h+1));   
% fase_med(2,:) = phiaf1(i:f,h*2+1);
% fase_med(3,:) = phiaf2(i:f,h*2+1);
% fase_med(4,:) = phiaf3(i:f,h*2+1);
% 
% 
% % atualiza amplitudes de referência para o calculo de erro das harmonicas
% ref_segH = ref_seg;
% ref_segH(2:4,i:f) = ref_seg(2:4,i:f)*amps(1,h+1);
% ref_segH(5:7,i:f) = (h+1)*ref_segH(5:7,i:f);
% 
% 
% [tve, amp_error, phase_error, freq_error] = calcula_erros(amp_med, fase_med, freq_final_fft, ref_segH(:,i:f));
% 
% % figure(3)
% % subplot(K,1,h+1), hold on, grid on%
% % plot(t_seg(i:f,1), tve(1,:), '*-')
% % plot(t_seg(i:f,1), tve(2,:), 'd-')
% % plot(t_seg(i:f,1), tve(3,:), 'x-')
% % plot(t_seg(i:f,1), tve(4,:), 's-')
% % %plot(t_seg(i:f,1), tve(5,:), '+-')
% % ylabel('TVE (%)'), xlabel('tempo (s)')
% 
% figure(3)
% subplot(K,1,h+1), hold on, grid on
% plot(t_seg(i:f,1), amp_med(1,:), '*-')
% plot(t_seg(i:f,1), amp_med(2,:), 'x-')
% plot(t_seg(i:f,1), amp_med(3,:), 'o-')
% plot(t_seg(i:f,1), amp_med(4,:), 's-')
% 
% 
% figure(4)
% subplot(K,1,h+1), hold on, grid on
% % plot(t_seg(i:f,1), amp_med(1,:), '*-')
% % plot(t_seg(i:f,1), amp_med(2,:), 'd-')
% % plot(t_seg(i:f,1), amp_med(3,:), 'x-')
% % plot(t_seg(i:f,1), amp_med(4,:), 's-')
% plot(t_seg(i:f,1), amp_error(1,:), '*-')
% plot(t_seg(i:f,1), amp_error(2,:), 'x-')
% plot(t_seg(i:f,1), amp_error(3,:), 'o-')
% plot(t_seg(i:f,1), amp_error(4,:), 's-')
% 
% 
% 
% 
% figure(5)
% subplot(K,1,h+1), hold on, grid on
% boxplot(amp_error(1:4,:)')
% 
% % figure(5)
% % subplot(K,1,h+1), hold on, grid on
% % plot(t_seg(i:f,1), fase_med(1,:), '*-')
% % plot(t_seg(i:f,1), fase_med(2,:), 'd-')
% % plot(t_seg(i:f,1), fase_med(3,:), 'x-')
% % plot(t_seg(i:f,1), fase_med(4,:), 's-')
% 
% end