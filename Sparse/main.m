%%% TESTES ESTIMADOR DFT-PRONY

clc
clear all
close all


set(0, 'defaultAxesFontSize',12);
set(0, 'defaultAxesFontName','times');




% ---------------------------------------------------------------------
% INICIALIZAÇÕES

N = 256;       % número de amostras
m = 256;       % passo de janelamento
Fs = 12800;    % Taxa de amostragem
f0 = 50;      % freq. fundamental (teórica)
noise = 0;    % nível de ruído em dB ( 0 = sem ruído)
T = 0.5;      % tempo total (em segundos)
param = 5;    % amort(w);    % freq. de modulação, por exemplo
type = 1;     % seleciona tipo de teste (ver gerador_sinais.m)

%freqs = [50];
freqs = [50 150 250 350 450];        % Vetor de frequências

%phases = zeros(3,2);
phases = [0 0 0 0 0;           % vetor de fases
         -2*pi/3 0 0 0 0;
         +2*pi/3 0 0 0 0];    
% phases(1,1) = deg2rad(0);
% phases(2,1) = deg2rad(-120);
% phases(3,1) = deg2rad(+120);
% % 
% phases(1,2) = deg2rad(0);
% phases(2,2) = deg2rad(-120);
% phases(3,2) = deg2rad(+120);
% 
% phases(1,3) = deg2rad(0);
% phases(2,3) = deg2rad(-120);
% phases(3,3) = deg2rad(+120);

amps = [1 0.7 0.5 0.3 0.1; 1 0.7 0.5 0.3 0.1; 1 0.7 0.5 0.3 0.1];
%  amps = [6 2.21 0.812 0.299;           % vetor de amplitudes
%          1 0.3 0 0;
%          1 0.3 0 0]; 

% ------------------------------------------------------------------------------------------------------------------
% GERAÇÃO DO SINAL

[Va, Vb, Vc, t, refs, mod] = gerador_sinais(f0, freqs, phases, amps,Fs, T,noise, type, param);

figure, hold on, grid on
% plot(t,mod*10)
plot(t,Va)
plot(t,Vb)
plot(t,Vc)
legend('Va','Vb','Vc'), title('Sinais')
xlabel('tempo (s)'), ylabel('amplitude')

%-------------------------------------------------------------------------------------------------------------------
% SEGMENTAÇÃO (JANELAMENTO)
[Va_seg, Vb_seg, Vc_seg, t_seg, ref_seg] = segmentador(Va, Vb, Vc, t, refs, N, m);



% % Plot amostragem (ciclo a ciclo)
% w = 4;   % numero de janelas no plot
% figure
% for i=1:w
%     subplot(w,1,i), hold on, grid on
%     plot(t_seg(i,:),Va_seg(i,:),'o-')
% end

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
% stem(gradeN,2*(Sa(1,:))/N)
% % stem(abs(Sa))
% % stem(abs(Sa))
% title('Amplitude (pu) - DFT')
% ylabel('pu'), xlabel('Frequência (Hz)')



%--------------------------------------------------------------------------------------------------------------------
% CÁLCULO DE COMPONENTES SIMÉTRICAS

bin = round(granN/f0)+1;       % local da componente fundamental no espectro

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


% % Plot fasores da componente fundamental
% figure
% subplot(2,3,1), hold on, grid on, plot(A0r, 's-'), plot(A0, 'o-'), title('Amplitudes (pu)'),ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,4), hold on, grid on, plot(rad2deg(phi0r), 's-'), plot(rad2deg(phi0), 'o-'), title('Fase (º)'), ylabel('º'), xlabel('Ciclo')
% subplot(2,3,2), hold on, grid on, plot(A1r, 's-'), plot(A1, 'o-'), title('Amplitudes (pu)'), ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,5), hold on, grid on, plot(rad2deg(phi1r), 's-'), plot(rad2deg(phi1), 'o-'), title('Fase (º)'), ylabel('º'), xlabel('Ciclo')
% subplot(2,3,3), hold on, grid on, plot(A2r, 's-'), plot(A2, 'o-'), title('Amplitudes (pu)'), ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,6), hold on, grid on, plot(rad2deg(phi2r), 's-'), plot(rad2deg(phi2), 'o-'), title('Fase (º)'), ylabel('º'), xlabel('Ciclo')


%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR DE FREQUÊNCIA - Derivada do ângulo do fasor de seq. positiva

% COM DADOS DE REFERÊNCIA (teste do estimador)
[freq_final_ref] = estimador_freq_deriv_ang(phi0r, phi1r, phi2r, f0, Fs, ref_seg, m, N);

% COM DADOS FFT
[freq_final_fft] = estimador_freq_deriv_ang(phi0, phi1, phi2, f0, Fs, ref_seg, m, N);

% figure, hold on, grid on
% plot(ref_seg(1,1:end), 'o-')
% plot(freq_final_ref, 'o-')
% plot(freq_final_fft, 'o-')
% title('Frequência (Hz)'), legend('Referência', 'Estimativa (ref)', 'Estimativa (DFT)');
% ylabel('Hz'), xlabel('Ciclo')

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR ESPARSO - DIC. FREQUÊNCIA - OMP_2

        % Fator de interpolação (expansão do dicionário)
        % Esparsidade

% Ajusta amplitude do espectro fourier
Sa = 2*Sa/N;
K = 1; 

P = 4;
H = N*P;
granH = (Fs/(H)) % Nova granularidade
gradeH = 0:granH:granH*((H)-1);

%load('C:\dicionarios_cs_dft\block_dic32u.mat')
%[Dic] = dic_freq(N, H);
load('Dic_A.mat')
%P = 1000;
omp = 1;
[Saf1, Aaf1, phiaf1, indx_af1] = estimador_freq_omp(Sa, K, P, omp, Fs, freq_final_fft, Dic);
%[Saf1, Aaf1, phiaf1, indx_af1] = estimador_freq_CORR(Sa, K, P, omp, Fs, freq_final_fft, Dic);
clear Dic % libera memoria


P = 10;
omp = 4;
[Saf2, Aaf2, phiaf2, indx_af2] = estimador_freq_omp_A(Sa, K, P, omp, Fs, f0, freq_final_fft);



%Aaf2 = Aaf2*2;
load('C:\dicionarios_cs_dft\block_dic32us.mat')
P = 10;
omp = 5;
[Saf3, Aaf3, phiaf3, indx_af3] = estimador_freq_omp_B(Sa, K, P, omp, Fs, f0, freq_final_fft, block);

clear block % libera memoria



%--------------------------------------------------------------------------------------------------------------------
% CÁLCULO DE ERROS

% índice da dft, alterar conforme incializações
hf(1) = 1 + 1; % fundamental

i=1;  % inicio
[f,~]=size(Aa); % fim (numero de ciclos)

for h=0:0;
% amplitudes
amp_med(1,:) = Aa(i:f,hf(h+1));    % dft
amp_med(2,:) = Aaf1(i:f,h*2+1);
% amp_med(3,:) = Aaf2(i:f,h*2+1); 
% amp_med(4,:) = Aaf3(i:f,h*2+1);

% fases
fase_med(1,:) = phia(i:f,hf(h+1));   
fase_med(2,:) = phiaf1(i:f,h*2+1);
% fase_med(3,:) = phiaf2(i:f,h*2+1);
% fase_med(4,:) = phiaf3(i:f,h*2+1);

[tve, amp_error, phase_error, freq_error] = calcula_erros(amp_med, fase_med, freq_final_fft, ref_seg(:,i:f));

end

%%%%%%%%%%%%%%%%%%%%%%%%% PLOT RESULTADOS

% % Plot fasores da componente fundamental
% figure
% subplot(2,1,1), hold on, grid on
% %plot(ref_seg(2,:),'o-')
% plot(t_seg(i:f,1), amp_med(1,:), '*-')
% plot(t_seg(i:f,1), amp_med(2,:), 'd-')
% plot(t_seg(i:f,1), amp_med(3,:), 'x-')
% plot(t_seg(i:f,1), amp_med(4,:), 's-')
% %plot(t_seg(i:f,1), amp_med(5,:), '+-')
% 
% ylabel('Amplitude (pu)'), xlabel('tempo (s)')
% subplot(2,1,2), hold on, grid on
% plot(t_seg(i:f,1),(fase_med(1,:)), '*-')
% plot(t_seg(i:f,1),(fase_med(2,:)), 'd-')
% plot(t_seg(i:f,1),(fase_med(3,:)), 'x-')
% plot(t_seg(i:f,1),(fase_med(4,:)), 's-')
% %plot(t_seg(i:f,1),(fase_med(5,:)), '+-')
% ylabel('Fase (rad)'), xlabel('tempo (s)')


figure, hold on, grid on
plot(t_seg(i:f,1), tve(1,:), '*-')
plot(t_seg(i:f,1), tve(2,:), 'x-')
% plot(t_seg(i:f,1), tve(3,:), 'o-')
% plot(t_seg(i:f,1), tve(4,:), 's-')
legend('DFT', 'Dic - Amp')
%legend('DFT', 'CS-DE', 'CS-DA', 'CS-DM')
ylabel('TVE (%)'), xlabel('tempo (s)')

% % Plot ERRO de fasores (DFT e Prony)
% figure  
% subplot(3,1,1), hold on, grid on%
% % for n=1:length(tve(:,1))
% %     plot(t_seg(i:f,1), tve(n,:), 'o-')
% % end
% plot(t_seg(i:f,1), tve(1,:), '*-')
% plot(t_seg(i:f,1), tve(2,:), 'x-')
% plot(t_seg(i:f,1), tve(3,:), 'o-')
% plot(t_seg(i:f,1), tve(4,:), 's-')
% legend('DFT', 'CS-Du', 'CS-Ddin', 'CS-Dmult')
% ylabel('TVE (%)'), xlabel('tempo (s)')
% 
% subplot(3,1,2), hold on, grid on
% % for n=1:length(tve(:,1))
% %     plot(t_seg(i:f,1), amp_error(n,:), 'o-')
% % end
% plot(t_seg(i:f,1),amp_error(1,:), '*-')
% plot(t_seg(i:f,1),amp_error(2,:), 'd-')
% plot(t_seg(i:f,1),amp_error(3,:), 'x-')
% plot(t_seg(i:f,1),amp_error(4,:), 's-')
% %plot(t_seg(i:f,1),amp_error(5,:), '+-')
% ylabel('Erro de amplitude (%)'), xlabel('tempo (s)')
% subplot(3,1,3), hold on, grid on
% % for n=1:length(tve(:,1))
% %     plot(t_seg(i:f,1), phase_error(n,:), 'o-')
% % end
% plot(t_seg(i:f,1),phase_error(1,:), '*-')
% plot(t_seg(i:f,1),phase_error(2,:), 'd-')
% plot(t_seg(i:f,1),phase_error(3,:), 'x-')
% plot(t_seg(i:f,1),phase_error(4,:), 's-')
% %plot(t_seg(i:f,1),phase_error(5,:), '+-')
% %legend('DFT', 'CS-Du', 'CS-Ddin', 'CS-Dmult')%, 'Freq. (OMP5)', 'Tempo (OMP1)', 'Tempo (OMP2)', 'Tempo (OMP5)')
% ylabel('Erro de fase (º)'), xlabel('tempo (s)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(tve(:,1));
   max_tve(i) = max(tve(i,:));
   max_amp_error(i) = max(amp_error(i,:)); 
   max_phi_error(i) = max(phase_error(i,:)); 
   
   std_tve(i) = std(tve(i,:));
   std_amp_error(i) = std(amp_error(i,:)); 
   std_phi_error(i) = std(phase_error(i,:)); 
end

% %%%%%%%%%%%%%%%%%%%%%%%%% ANÁLISE DAS HARMÔNICAS
% 
% % índice da dft, alterar conforme incializações
% hf(1) = 1 + 1; % fundamental
% hf(2) = 1 + 3; % 3th harmonica
% hf(3) = 1 + 5; % 5th harmonica
% 
% i=1;  % inicio
% [f,~]=size(Aa); % fim (numero de ciclos)
% 
% 
% for h=0:0
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
% 
% 
% [tve, amp_error, phase_error, freq_error] = calcula_erros(amp_med, fase_med, freq_final_fft, ref_segH(:,i:f));
% 
% figure(3)
% subplot(3,1,h+1), hold on, grid on%
% plot(t_seg(i:f,1), tve(1,:), '*-')
% plot(t_seg(i:f,1), tve(2,:), 'd-')
% plot(t_seg(i:f,1), tve(3,:), 'x-')
% plot(t_seg(i:f,1), tve(4,:), 's-')
% %plot(t_seg(i:f,1), tve(5,:), '+-')
% ylabel('TVE (%)'), xlabel('tempo (s)')
% 
% % if h==0
% %     figure
% % end
% figure(4)
% subplot(3,1,h+1), hold on, grid on
% plot(t_seg(i:f,1), amp_med(1,:), '*-')
% plot(t_seg(i:f,1), amp_med(2,:), 'd-')
% plot(t_seg(i:f,1), amp_med(3,:), 'x-')
% plot(t_seg(i:f,1), amp_med(4,:), 's-')
% 
% end