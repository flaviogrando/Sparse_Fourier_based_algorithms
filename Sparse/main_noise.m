%%% TESTES ESTIMADOR CS-DFT
%
% 

clc
clear all
close all

load('E:\block_dic256u.mat')
load('E:\block_dic256us.mat')

set(0, 'defaultAxesFontSize',12);
set(0, 'defaultAxesFontName','times');
wb = waitbar(0,'Progresso');

ff = 30:5:70; % cria vetor para deslocar o ponto da grade

for w=1:length(ff)

clc
waitbar(w/length(ff));
% ---------------------------------------------------------------------
% INICIALIZAÇÕES

N = 256;      % número de amostras
m = 256;      % passo de janelamento
Fs = 12800;    % Taxa de amostragem
f0 = 50;      % freq. fundamental (teórica)
noise = ff(w)   % nível de ruído em dB ( 0 = sem ruído)
T = 1;        % tempo total (em segundos)
param = 5;    % amort(w);    % freq. de modulação, por exemplo
type = 0;     % seleciona tipo de teste (ver gerador_sinais.m)

freqs = [50.025]; % Vetor de frequências

%phases = zeros(3,2);   % vetor de fases
phases(1,1) = deg2rad(0); 
phases(2,1) = deg2rad(-120);
phases(3,1) = deg2rad(+120);

amps = [1; 1; 1;]; % vetor de amplitudes

% freqs = [49.5 99 148.51 198.01]; 
% amps = [6 2.21 0.812 0.299;           
%         1 0.3 0 0;
%         1 0.3 0 0]; 
% phases = [3.96 4.85 4.75 3.11;         
%          -2*pi/3 0 0 0;
%          +2*pi/3 0 0 0];  

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
% stem(2*abs(Sa(1,:))/N)
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


P = 1000;    % Fator de interpolação (expansão do dicionário)
K = 1;         % Esparsidade
H = P*N;
granH = (Fs/(H)); % Nova granularidade
gradeH = 0:granH:granH*((H)-1);
% Ajusta amplitude do espectro fourier
Sa = 2*Sa/N;


%[Dic] = dic_freq(N, H);
omp = 4;
 
[Saf1, Aaf1, phiaf1, indx_af1] = estimador_freq_omp(Sa, K, P, omp, Fs, freq_final_fft, Dic);


P = 10;
omp = 4;
[Saf2, Aaf2, phiaf2, indx_af2] = estimador_freq_omp_A(Sa, K, P, omp, Fs, f0, freq_final_fft);



P = 10;
omp = 5;
[Saf3, Aaf3, phiaf3, indx_af3] = estimador_freq_omp_B(Sa, K, P, omp, Fs, f0, freq_final_fft, block);

%clear block % libera memoria


%--------------------------------------------------------------------------------------------------------------------
% CÁLCULO DE ERROS

% define indices (recorte dos dados)
%[num_win ] = size(Aa);
i=1;  % inicio
%[f,~] = size(Aa); % fim (numero de ciclos)
f=49;

indices_min(1,w) = 0; 
indices_min(2,w) = min(indx_af1(1,:)); 
indices_min(3,w) = min(indx_af2(1,:)); 
indices_min(4,w) = min(indx_af3(1,:)); 

indices_max(1,w) = 0; 
indices_max(2,w) = max(indx_af1(1,:)); 
indices_max(3,w) = max(indx_af2(1,:)); 
indices_max(4,w) = max(indx_af3(1,:)); 


% amplitudes
amp_med(1,:) = Aa(i:f,2);    % dft fundamental
amp_med(2,:) = Aaf1(i:f,1);
amp_med(3,:) = Aaf2(i:f,1);  
amp_med(4,:) = Aaf3(i:f,1);
%amp_med(5,w) = Aaf4(i:f,1);
 
% fases
fase_med(1,:) = phia(i:f,2);   
fase_med(2,:) = phiaf1(i:f,1);
fase_med(3,:) = phiaf2(i:f,1);
fase_med(4,:) = phiaf3(i:f,1);
%fase_med(5,w) = phiaf4(i:f,1);

% frequência
%freq_med(1,:) = ref_seg(1,:);
freq_med(1,:) = freq_final_ref(1,i:f);
freq_med(2,:) = freq_final_fft(1,i:f);
% frequência


% % novo conjunto de referÊncias segmentadas
 %new_ref_seg(:,w) = ref_seg(:,i:f);


[tve, amp_error, phase_error, freq_error] = calcula_erros(amp_med, fase_med, freq_med, ref_seg(:,:));


    % ERROS EM FUNÇÃO DO RUÍDO
    for i=1:length(tve(1,:)) 
        tve_dft(:,w) = tve(1,:);        
        tve_cs_de(:,w) = tve(2,:); 
        tve_cs_da(:,w) = tve(3,:);     
        tve_cs_dm(:,w) = tve(4,:);
    end

end
close(wb)

figure
subplot(4,1,1), hold on, grid on%
boxplot(tve_dft, ff)
xlabel('SNR (dB)'), ylabel('TVE (%)'), legend('DFT')
subplot(4,1,2), hold on, grid on%
boxplot(tve_cs_de, ff)
xlabel('SNR (dB)'), ylabel('TVE (%)'), legend('CS-DE')
subplot(4,1,3), hold on, grid on%
boxplot(tve_cs_da, ff)
xlabel('SNR (dB)'), ylabel('TVE (%)'), legend('CS-DA')
subplot(4,1,4), hold on, grid on%
boxplot(tve_cs_dm, ff)
xlabel('SNR (dB)'), ylabel('TVE (%)'), legend('CS-DM')