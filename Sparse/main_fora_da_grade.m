%%% TESTES ESTIMADOR CS-DFT
%
% 

clc
%clear all
%close all


set(0, 'defaultAxesFontSize',12);
set(0, 'defaultAxesFontName','times');
wb = waitbar(0,'Progresso');

ff=-0.15:0.005:0.15; % cria vetor para deslocar o ponto da grade

for w=1:length(ff)
clc
waitbar(w/length(ff));
% ---------------------------------------------------------------------
% INICIALIZAÇÕES

N = 256;      % número de amostras
m = 1;      % passo de janelamento
Fs = 12800;    % Taxa de amostragem.
f0 = 60;      % freq. fundamental (teórica)
noise = 50;   % nível de ruído em dB ( 0 = sem ruído)
T = 1*N/Fs + 1/Fs;        % tempo total (em segundos)
param = 5;    % amort(w);    % freq. de modulação, por exemplo
type = 0;     % seleciona tipo de teste (ver gerador_sinais.m)

freqs = [50 + ff(w)]
% freqs = [49.5 99 148.51 198.01];        % Vetor de frequências

%phases = zeros(3,2);
% phases = [3.96 4.85 4.75 3.11;           % vetor de fases
%          -2*pi/3 0 0 0;
%          +2*pi/3 0 0 0];    
phases(1,1) = deg2rad(0);
phases(2,1) = deg2rad(-120);
phases(3,1) = deg2rad(+120);

amps = [1; 1; 1;];
% amps = [6 2.21 0.812 0.299;           % vetor de amplitudes
%         1 0.3 0 0;
%         1 0.3 0 0]; 

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

load('E:\block_dic256u.mat')
%[Dic] = dic_freq(N, H);
omp = 4;
 
[Saf1, Aaf1, phiaf1, indx_af1] = estimador_freq_omp(Sa, K, P, omp, Fs, freq_final_fft, Dic);


P = 10;
omp = 4;
[Saf2, Aaf2, phiaf2, indx_af2] = estimador_freq_omp_A(Sa, K, P, omp, Fs, f0, freq_final_fft);


load('E:\block_dic256us.mat')
P = 10;
omp = 5;
[Saf3, Aaf3, phiaf3, indx_af3] = estimador_freq_omp_B(Sa, K, P, omp, Fs, f0, freq_final_fft, block);

%clear block % libera memoria





%--------------------------------------------------------------------------------------------------------------------
% CÁLCULO DE ERROS



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
amp_med(1,w) = Aa(i:f,2);    % dft fundamental
amp_med(2,w) = Aaf1(i:f,1);
amp_med(3,w) = Aaf2(i:f,1);  
amp_med(4,w) = Aaf3(i:f,1);
%amp_med(5,w) = Aaf4(i:f,1);
 
% fases
fase_med(1,w) = phia(i:f,2);   
fase_med(2,w) = phiaf1(i:f,1);
fase_med(3,w) = phiaf2(i:f,1);
fase_med(4,w) = phiaf3(i:f,1);
%fase_med(5,w) = phiaf4(i:f,1);

% frequência
%freq_med(1,:) = ref_seg(1,:);
freq_med(1,w) = freq_final_ref(1,i:f);
freq_med(2,w) = freq_final_fft(1,i:f);


% novo conjunto de referÊncias segmentadas
 new_ref_seg(:,w) = ref_seg(:,i:f);
 end
close(wb);

 [tve, amp_error, phase_error, freq_error] = calcula_erros(amp_med, fase_med, freq_med, new_ref_seg(:,:));



% figure
% subplot(2,1,1), hold on, grid on
% %plot(ff, indices_min(1,:), '*-')
% plot(ff, indices_min(2,:), 'd-')
% plot(ff, indices_min(3,:), 'x-')
% plot(ff, indices_min(4,:), 's-')
% title('indices')
% subplot(2,1,2), hold on, grid on
% %plot(ff, indices_max(1,:), '*-')
% plot(ff, indices_max(2,:), 'd-')
% plot(ff, indices_max(3,:), 'x-')
% plot(ff, indices_max(4,:), 's-')

% % Plot fasores da componente fundamental
% figure
% subplot(2,1,1), hold on, grid on
% %plot(ref_seg(2,:),'o-')
% plot(ff, amp_med(1,:), '*-')
% plot(ff, amp_med(2,:), 'x-')
% plot(ff, amp_med(3,:), 'o-')
% plot(ff, amp_med(4,:), 's-')
% %plot(ff, amp_med(5,:), '+-')
% ylabel('Amplitude (pu)'), xlabel('\Delta f')
% subplot(2,1,2), hold on, grid on
% plot(ff,(fase_med(1,:)), '*-')
% plot(ff,(fase_med(2,:)), 'x-')
% plot(ff,(fase_med(3,:)), 'o-')
% plot(ff,(fase_med(4,:)), 's-')
% %plot(ff,(fase_med(5,:)), '+-')
% ylabel('Fase (rad)'), xlabel('\Delta f')
% legend('DFT', 'CS')
% %legend('DFT', 'Orig', 'Desl', 'Desl', 'Desl. otm');

figure  
hold on, grid on
plot(ff, tve(1,:), '*-')
plot(ff, tve(2,:), 'x-')
plot(ff, tve(3,:), 'o-')
plot(ff, tve(4,:), 's-')
ylabel('TVE (%)'), xlabel('\Deltaf')
legend('DFT', 'CS-DE', 'CS-DA', 'CS-DM')


% figure  
% hold on, grid on
% plot(ff, tve(1,:), '*-')
% plot(ff, tve(2,:))
% plot(ff, tve(3,:))
% plot(ff, tve(4,:))
% %plot(ff, tve(5,:))
% ylabel('TVE (%)'), xlabel('\Deltaf1')
% legend('DE (fe)', 'DE (OMP)', 'DA (OMP)', 'DM (OMP)')
% %legend('DE (fe)', 'DE (OMP)', 'DA (OMP)', 'DM (OMP)')


% % Plot ERRO de fasores (DFT e Prony)
% figure  
% subplot(3,1,1), hold on, grid on%
% % for n=1:length(tve(:,1))
% %     plot(ff, tve(n,:), 'o-')
% % end
% plot(ff, tve(1,:), '*-')
% plot(ff, tve(2,:), 'x-')
% plot(ff, tve(3,:), 'o-')
% plot(ff, tve(4,:), 's-')
% %plot(ff, tve(5,:), '+-')
% ylabel('TVE (%)'), xlabel('\Deltaf')
% subplot(3,1,2), hold on, grid on
% % for n=1:length(tve(:,1))
% %     plot(ff, amp_error(n,:), 'o-')
% % end
% plot(ff,amp_error(1,:), '*-')
% plot(ff,amp_error(2,:), 'x-')
% plot(ff,amp_error(3,:), 'o-')
% plot(ff,amp_error(4,:), 's-')
% %plot(ff,amp_error(5,:), '+-')
% ylabel('Erro de amplitude (%)'), xlabel('\Deltaf')
% subplot(3,1,3), hold on, grid on
% % for n=1:length(tve(:,1))
% %     plot(ff, phase_error(n,:), 'o-')
% % end
% plot(ff,phase_error(1,:), '*-')
% plot(ff,phase_error(2,:), 'x-')
% plot(ff,phase_error(3,:), 'o-')
% plot(ff,phase_error(4,:), 's-')
% %plot(ff,phase_error(5,:), '+-')
% %legend('DFT', 'Freq. (OMP1)', 'Freq. (OMP2)', 'Freq. (OMP5)', 'Tempo (OMP1)', 'Tempo (OMP2)', 'Tempo (OMP5)')
% ylabel('Erro de fase (º)'), xlabel('\Deltaf')
% legend('DFT', 'CS-D_{u}', 'CS-D_{din}', 'CS-D_{mult}')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:4;
   max_tve(i) = max(tve(i,:));
   max_amp_error(i) = max(amp_error(i,:)); 
   max_phi_error(i) = max(phase_error(i,:)); 
end
%