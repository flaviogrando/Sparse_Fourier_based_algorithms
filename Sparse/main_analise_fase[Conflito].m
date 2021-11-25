%%% TESTES ESTIMADOR CS-DFT
%
% 

clc
clear all
close all


set(0, 'defaultAxesFontSize',12);
set(0, 'defaultAxesFontName','times');
wb = waitbar(0,'Progresso');

ff=-175:5:180; % cria vetor para deslocar o ponto da grade

for w=1:length(ff)


    waitbar(w/length(ff));
    
% ---------------------------------------------------------------------
% INICIALIZA��ES

N = 128;      % n�mero de amostras
m = 1;      % passo de janelamento
Fs = 6400;    % Taxa de amostragem
f0 = 50;      % freq. fundamental (te�rica)
noise = 00;   % n�vel de ru�do em dB ( 0 = sem ru�do)
T = 0.02+1/Fs;        % tempo total (em segundos)
param = 5;    % amort(w);    % freq. de modula��o, por exemplo
type = 0;     % seleciona tipo de teste (ver gerador_sinais.m)

freqs = [52.5];
% freqs = [49.5 99 148.51 198.01];        % Vetor de frequ�ncias

%phases = zeros(3,2);
% phases = [3.96 4.85 4.75 3.11;           % vetor de fases
%          -2*pi/3 0 0 0;
%          +2*pi/3 0 0 0];    
phases(1,1) = deg2rad(0) + deg2rad(ff(w));
phases(2,1) = deg2rad(-120) + deg2rad(ff(w));
phases(3,1) = deg2rad(+120) + deg2rad(ff(w));

amps = [1; 1; 1;];
% amps = [6 2.21 0.812 0.299;           % vetor de amplitudes
%         1 0.3 0 0;
%         1 0.3 0 0]; 

% ------------------------------------------------------------------------------------------------------------------
% GERA��O DO SINAL

[Va, Vb, Vc, t, refs, mod] = gerador_sinais(f0, freqs, phases, amps,Fs, T,noise, type, param);

% figure, hold on, grid on
% % plot(t,mod*10)
% plot(t,Va)
% plot(t,Vb)
% plot(t,Vc)
% legend('Va','Vb','Vc'), title('Sinais')
% xlabel('tempo (s)'), ylabel('amplitude')

%-------------------------------------------------------------------------------------------------------------------
% SEGMENTA��O (JANELAMENTO)
[Va_seg, Vb_seg, Vc_seg, t_seg, ref_seg] = segmentador(Va, Vb, Vc, t, refs, N, m);

% % Plot amostragem (ciclo a ciclo)
% w = 4;   % numero de janelas no plot
% figure
% for i=1:w
%     subplot(w,1,i), hold on, grid on
%     plot(t_seg(i,:),Va_seg(i,:),'o-')
% end

% figure %Plot refer�ncias
% subplot(3,1,1), hold on, grid on
% plot(t_seg(:,1), ref_seg(1,:), 'o-')
% xlabel('tempo (s)'), ylabel('freq (Hz)'), title('Refer�ncias (segmentada)')
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

% granN = (Fs/N);
% eixo_x = 0:granN:granN*(N-1);
% % 
% % Plot espectro
% figure
% hold on, grid on
% stem(eixo_x,2*abs(Sa(1,:))/N)
% % stem(abs(Sa))
% % stem(abs(Sa))
% title('Amplitude (pu) - DFT')
% ylabel('pu'), xlabel('Frequ�ncia (Hz)')



%--------------------------------------------------------------------------------------------------------------------
% C�LCULO DE COMPONENTES SIM�TRICAS

% COM DADOS FFT
[A0, A1, A2, phi0, phi1, phi2] = comp_simetricas(Aa, Ab, Ac, phia, phib, phic);

% COM DADOS DE REFER�NCIA (para teste do estimador)
[A0r, A1r, A2r, phi0r, phi1r, phi2r] = comp_simetricas(ref_seg(2,:), ref_seg(3,:), ref_seg(4,:),ref_seg(5,:), ref_seg(6,:), ref_seg(7,:));


% % Plot fasores da componente fundamental
% figure
% subplot(2,3,1), hold on, grid on, plot(A0r, 's-'), plot(A0, 'o-'), title('Amplitudes (pu)'),ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,4), hold on, grid on, plot(rad2deg(phi0r), 's-'), plot(rad2deg(phi0), 'o-'), title('Fase (�)'), ylabel('�'), xlabel('Ciclo')
% subplot(2,3,2), hold on, grid on, plot(A1r, 's-'), plot(A1, 'o-'), title('Amplitudes (pu)'), ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,5), hold on, grid on, plot(rad2deg(phi1r), 's-'), plot(rad2deg(phi1), 'o-'), title('Fase (�)'), ylabel('�'), xlabel('Ciclo')
% subplot(2,3,3), hold on, grid on, plot(A2r, 's-'), plot(A2, 'o-'), title('Amplitudes (pu)'), ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,6), hold on, grid on, plot(rad2deg(phi2r), 's-'), plot(rad2deg(phi2), 'o-'), title('Fase (�)'), ylabel('�'), xlabel('Ciclo')

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR DE FREQU�NCIA - Derivada do �ngulo do fasor de seq. positiva

% COM DADOS DE REFER�NCIA (teste do estimador)
[freq_final_ref] = estimador_freq_deriv_ang(phi0r, phi1r, phi2r, f0, Fs, ref_seg, m, N);

% COM DADOS FFT
[freq_final_fft] = estimador_freq_deriv_ang(phi0, phi1, phi2, f0, Fs, ref_seg, m, N);

% figure, hold on, grid on
% plot(ref_seg(1,1:end), 'o-')
% plot(freq_final_ref, 'o-')
% plot(freq_final_fft, 'o-')
% title('Frequ�ncia (Hz)'), legend('Refer�ncia', 'Estimativa (ref)', 'Estimativa (DFT)');
% ylabel('Hz'), xlabel('Ciclo')

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR ESPARSO - DIC. FREQU�NCIA - OMP_2

P = 20;        % Fator de interpola��o (expans�o do dicion�rio)
K = 2;         % Esparsidade
granH = (Fs/(P*N)); % Nova granularidade
% Ajusta amplitude do espectro fourier
Sa = 2*Sa/N;



omp = 1;
[Saf1, Aaf1, phiaf1, indx_af1] = estimador_freq_omp(Sa, K, P, omp);

omp = 2;
[Saf2, Aaf2, phiaf2, indx_af2] = estimador_freq_omp(Sa, K, P, omp);

omp = 3;
[Saf3, Aaf3, phiaf3, indx_af3] = estimador_freq_omp(Sa, K, P, omp);


% omp = 1;
% [Saf4, Aaf4, phiaf4, indx_af4] = estimador_freq_omp_A(Sa, K, P, omp, Fs, f0, freq_final_fft);
% 
% omp = 2;
% [Saf5, Aaf5, phiaf5, indx_af5] = estimador_freq_omp_A(Sa, K, P, omp, Fs, f0, freq_final_fft);
% 
% omp = 3;
% [Saf6, Aaf6, phiaf6, indx_af6] = estimador_freq_omp_A(Sa, K, P, omp, Fs, f0, freq_final_fft);


% omp = 1;
% [Saf1, Aaf1, phiaf1, indx_af1] = estimador_time_omp(Sa, K, P, omp);
% 
% omp = 2;
% [Saf2, Aaf2, phiaf2, indx_af2] = estimador_time_omp(Sa, K, P, omp);
% 
% omp = 3;
% [Saf3, Aaf3, phiaf3, indx_af3] = estimador_time_omp(Sa, K, P, omp);


% omp = 1;
% [Saf4, Aaf4, phiaf4, indx_af4] = estimador_time_omp_A(Sa, K, P, omp, Fs, f0, freq_final_fft);
% 
% omp = 2;
% [Saf5, Aaf5, phiaf5, indx_af5] = estimador_time_omp_A(Sa, K, P, omp, Fs, f0, freq_final_fft);
% 
% omp = 3;
% [Saf6, Aaf6, phiaf6, indx_af6] = estimador_time_omp_A(Sa, K, P, omp, Fs, f0, freq_final_fft);



%--------------------------------------------------------------------------------------------------------------------
% C�LCULO DE ERROS



% define indices (recorte dos dados)
%[num_win ] = size(Aa);
i=1;  % inicio
f=1; %length(Aa); % fim (numero de ciclos)

indices_min(1,w) = 0; 
indices_min(2,w) = min(indx_af1(1,:)); 
indices_min(3,w) = min(indx_af2(1,:)); 
indices_min(4,w) = min(indx_af3(1,:)); 

indices_max(1,w) = 0; 
indices_max(2,w) = max(indx_af1(1,:)); 
indices_max(3,w) = max(indx_af2(1,:)); 
indices_max(4,w) = max(indx_af3(1,:)); 


% amplitudes
amp_med(1,w) = Aa(1,i:f);    % dft
amp_med(2,w) = Aaf1(1,i:f);
amp_med(3,w) = Aaf2(1,i:f);  
amp_med(4,w) = Aaf3(1,i:f);
% amp_med(5,w) = Aaf7(1,i:f);
% amp_med(6,w) = Aat2(1,i:f);  
% amp_med(7,w) = Aat5(1,i:f);  
% fases
fase_med(1,w) = phia(1,i:f);   
fase_med(2,w) = phiaf1(1,i:f);
fase_med(3,w) = phiaf2(1,i:f);
fase_med(4,w) = phiaf3(1,i:f);
% fase_med(5,w) = phiaf7(1,i:f);
% fase_med(6,w) = phiat2(1,i:f);
% fase_med(7,w) = phiat5(1,i:f);
% frequ�ncia
% %freq_med(1,:) = ref_seg(1,:);
freq_med(1,w) = freq_final_ref(1,i:f);
freq_med(2,w) = freq_final_fft(1,i:f);


% novo conjunto de refer�ncias segmentadas
 new_ref_seg(:,w) = ref_seg(:,i:f);
 end
close(wb);

 [tve, amp_error, phase_error, freq_error] = calcula_erros(amp_med, fase_med, freq_med, new_ref_seg(:,:));



% figure
% subplot(2,1,1), hold on, grid on
% plot(ff, indices_min(1,:), '*-')
% plot(ff, indices_min(2,:), 'd-')
% plot(ff, indices_min(3,:), 'x-')
% plot(ff, indices_min(4,:), 's-')
% title('indices')
% subplot(2,1,2), hold on, grid on
% plot(ff, indices_max(1,:), '*-')
% plot(ff, indices_max(2,:), 'd-')
% plot(ff, indices_max(3,:), 'x-')
% plot(ff, indices_max(4,:), 's-')

% Plot fasores da componente fundamental
figure
subplot(2,1,1), hold on, grid on
%plot(ref_seg(2,:),'o-')
%plot(ff, amp_med(1,:), '*-')
plot(ff, amp_med(2,:), 'd-')
plot(ff, amp_med(3,:), 'x-')
plot(ff, amp_med(4,:), 's-')
%plot(ff, amp_med(5,:), '+-')
ylabel('Amplitude (pu)'), xlabel('\Delta f')
subplot(2,1,2), hold on, grid on
%plot(ff,(fase_med(1,:)), '*-')
plot(ff,(fase_med(2,:)), 'd-')
plot(ff,(fase_med(3,:)), 'x-')
plot(ff,(fase_med(4,:)), 's-')
%plot(ff,(fase_med(5,:)), '+-')
ylabel('Fase (rad)'), xlabel('degrees')
legend('DFT', 'Orig', 'Desl', 'Desl', 'Desl. otm');




% Plot ERRO de fasores (DFT e Prony)
figure  
subplot(3,1,1), hold on, grid on%
% for n=1:length(tve(:,1))
%     plot(ff, tve(n,:), 'o-')
% end
%plot(ff, tve(1,:), '*-')
plot(ff, tve(2,:), 'd-')
plot(ff, tve(3,:), 'x-')
plot(ff, tve(4,:), 's-')
%plot(ff, tve(5,:), '+-')

ylabel('TVE (%)'), xlabel('degrees')
subplot(3,1,2), hold on, grid on
% for n=1:length(tve(:,1))
%     plot(ff, amp_error(n,:), 'o-')
% end
%plot(ff,amp_error(1,:), '*-')
plot(ff,amp_error(2,:), 'd-')
plot(ff,amp_error(3,:), 'x-')
plot(ff,amp_error(4,:), 's-')
%plot(ff,amp_error(5,:), '+-')


ylabel('Erro de amplitude (%)'), xlabel('degrees')
subplot(3,1,3), hold on, grid on
% for n=1:length(tve(:,1))
%     plot(ff, phase_error(n,:), 'o-')
% end
%plot(ff,phase_error(1,:), '*-')
plot(ff,phase_error(2,:), 'd-')
plot(ff,phase_error(3,:), 'x-')
plot(ff,phase_error(4,:), 's-')
%plot(ff,phase_error(5,:), '+-')

%legend('DFT', 'Freq. (OMP1)', 'Freq. (OMP2)', 'Freq. (OMP5)', 'Tempo (OMP1)', 'Tempo (OMP2)', 'Tempo (OMP5)')
ylabel('Erro de fase (�)'), xlabel('degrees')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:2;
   max_tve(i) = max(tve(i,:));
   max_amp_error(i) = max(amp_error(i,:)); 
   max_phi_error(i) = max(phase_error(i,:)); 
end
%