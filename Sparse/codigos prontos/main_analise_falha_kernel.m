%%% TESTES ESTIMADOR ESPARSO

clc
clear all
close all


set(0, 'defaultAxesFontSize',12);
set(0, 'defaultAxesFontName','times');

%Ddk = zeros(100,2000);
%fail = zeros(1,20);
for i=1:20
    fail(i) = 10^(-i);
end


 for w=1:length(fail)


%amort_vetor = ones(50)*amort(w);   % Vetor de referência para amortecimento

% ---------------------------------------------------------------------
% INICIALIZAÇÕES

N = 160;      % número de amostras
m = 200;      % passo de janelamento
Fs = 4000;   % Taxa de amostragem
f0 = 50;      % freq. fundamental (teórica)
noise = 0;    % nível de ruído em dB ( 0 = sem ruído)
T = N/Fs;        % tempo total (em segundos)
param = 0;   % amort(w);    % freq. de modulação, por exemplo
type = 0;     % seleciona tipo de teste (ver gerador_sinais.m)

freqs = [50];        % Vetor de frequências

phases = zeros(3,2);
% phases = [0 0;           % vetor de fases
%          -2*pi/3 0;
%          +2*pi/3 0];    
phases(1,1) = deg2rad(0);
phases(2,1) = deg2rad(-120);
phases(3,1) = deg2rad(+120);

amps = [1 0.3;           % vetor de amplitudes
        1 0.3;
        1 0.3]; 

% ------------------------------------------------------------------------------------------------------------------
% GERAÇÃO DO SINAL
[Va, Vb, Vc, t, refs, mod] = gerador_sinais(f0, freqs, phases, amps,Fs, T,noise, type, param);

% figure, hold on, grid on
% %plot(t,mod*10)
% plot(t,Va)
% plot(t,Vb)
% plot(t,Vc)
% legend('Va','Vb','Vc'), title('Sinais')
% xlabel('tempo (s)'), ylabel('amplitude')

%-------------------------------------------------------------------------------------------------------------------
% SEGMENTAÇÃO (JANELAMENTO)
[Va_seg, Vb_seg, Vc_seg, t_seg, ref_seg] = segmentador(Va, Vb, Vc, t, refs, N, m);

% % Plot amostragem (ciclo a ciclo)
% w = 1;   % numero de janelas no plot
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

gran1 = (Fs/N);
eixo_x = 0:gran1:gran1*(N-1);

% % Plot espectro
% figure
% hold on, grid on
% stem(eixo_x,2*abs(Sa(1,:))/N)
% % stem(abs(Sa))
% % stem(abs(Sa))
% title('Amplitude (pu) - DFT')
% ylabel('pu'), xlabel('Frequência (Hz)')



%--------------------------------------------------------------------------------------------------------------------
% CÁLCULO DE COMPONENTES SIMÉTRICAS

% COM DADOS FFT
[A0, A1, A2, phi0, phi1, phi2] = comp_simetricas(Aa, Ab, Ac, phia, phib, phic);

% COM DADOS DE REFERÊNCIA (para teste do estimador)
[A0r, A1r, A2r, phi0r, phi1r, phi2r] = comp_simetricas(ref_seg(2,:), ref_seg(3,:), ref_seg(4,:),ref_seg(5,:), ref_seg(6,:), ref_seg(7,:));


% % Plot fasores da componente fundamental
% figure
% 
% subplot(2,3,1), hold on, grid on
% plot(A0r, 's-')
% plot(A0, 'o-')
% title('Amplitudes (pu)')
% ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,4), hold on, grid on
% plot(rad2deg(phi0r), 's-')
% plot(rad2deg(phi0), 'o-')
% 
% title('Fase (º)')
% ylabel('º'), xlabel('Ciclo')
% 
% subplot(2,3,2), hold on, grid on
% plot(A1r, 's-')
% plot(A1, 'o-')
% title('Amplitudes (pu)')
% ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,5), hold on, grid on
% plot(rad2deg(phi1r), 's-')
% plot(rad2deg(phi1), 'o-')
% title('Fase (º)')
% ylabel('º'), xlabel('Ciclo')
% 
% subplot(2,3,3), hold on, grid on
% plot(A2r, 's-')
% plot(A2, 'o-')
% title('Amplitudes (pu)')
% ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,6), hold on, grid on
% plot(rad2deg(phi2r), 's-')
% plot(rad2deg(phi2), 'o-')
% title('Fase (º)')
% ylabel('º'), xlabel('Ciclo')

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
% DIRICHLET KERNEL - expanded

P = 10; % fator de interpolação
H = N*P;


%Ddk = zeros(N,H);  % Dirichlet kernel
for l=0:H-1      % colunas
    for k=0:N-1  % linhas             
        arg = k/N - ((l+fail(w))/H); % argumento do kernel 
        %arg = k/N - ((l)/H); % argumento do kernel 
        Ddk(k+1,l+1) = (sin(pi*N*(arg))/(N*sin(pi*(arg))))*exp(-1i*pi*(N-1)*(arg)); % eq. (3)-1
        %Ddk(k+1,l+1) = (sin(pi*N*(arg))/(N*sin(pi*(arg)))+fail(w))*exp(-1i*pi*(N-1)*(arg));
    end
end

% figure
% mesh(abs(Ddk))
% xlabel('columns');
% ylabel('rows');
% title('Dirichlet kernel');

%----------------------------------------------------------------------------------------------------

Sa = 2*Sa/N;

% Algoritmo 1
m = 2;
errorGoal =  1e-10;
tic;
[spectro1, indice1]=OMP_phasor(Ddk, Sa',m,0,errorGoal);
toc;

% Algoritmo 2
tic;
[spectro2] = OMP_2(m,Sa,Ddk);
toc;

[spectro3] = omp_3(Ddk,Sa',m);

gran2 = (Fs/H);
eixo_x = 0:gran2:gran2*(H-1);

f1 = indice1(1,1); % indice da freq fundamental

Asa = abs(spectro1(f1,1));
phisa = angle(spectro1(f1,1));

Asa2 = abs(spectro2(1,f1));
phisa2 = angle(spectro2(1,f1));


% figure, hold on, grid on
% stem(eixo_x, Asa)

% Seleção automática da componente fundamental
[num_wins,win_size] = size(Aa);
gran = Fs/win_size/f0;       % granularidade
bin = round(1/gran)+1;       % local da componente fundamental no espectro



%--------------------------------------------------------------------------------------------------------------------
% CÁLCULO DE ERROS

% % Seleção automática da componente fundamental
% [num_wins,win_size] = size(Aa);
% gran = Fs/win_size/f0;       % granularidade
% bin = round(1/gran)+1;       % local da componente fundamental no espectro

% define indices (recorte dos dados)

i=1;  % inicio
f=1;%floor(T/(1/f0))-1; % fim (numero de ciclos)

% amplitudes
amp_med(1,w) = Aa(1,i:f);    % dft
amp_med(2,w) = Asa(1,i:f);  % com freq de ref
amp_med(3,w) = Asa2(1,i:f);  % com freq de ref
% fases
fase_med(1,w) = phia(1,i:f);   
fase_med(2,w) = phisa(1,i:f);
fase_med(3,w) = phisa2(1,i:f);

% frequência
%freq_med(1,:) = ref_seg(1,:);
freq_med(1,w) = freq_final_ref(1,i:f);
freq_med(2,w) = freq_final_fft(1,i:f);


% novo conjunto de referÊncias segmentadas
 new_ref_seg(:,w) = ref_seg(:,i:f);
 end


 [tve, amp_error, phase_error, freq_error] = calcula_erros(amp_med, fase_med, freq_med, new_ref_seg(:,:));
 
 
 
% % Plot erro frequência
% figure
% subplot(2,1,1),hold on, grid on  
% plot(ref_seg(1,:), 'o-')
% plot(freq_final, 'o-')
% plot(freq_final_prony, 'o-')
% plot(freq_final_ref, '*-')
% title('Frequência (Hz)'), legend('Ref', 'DFT', 'Prony');
% ylabel('Hz'), xlabel('Ciclo')
% subplot(2,1,2),hold on, grid on  
% plot(freq_error(1,:), 'o-')
% plot(freq_error(2,:), 'o-')
% plot(freq_error(3,:), 'o-')
% title('Erro de Frequência (Hz)'), legend('DFT', 'Prony');
% ylabel('Hz'), xlabel('Ciclo')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gambiarra pro t_seg
i=1;  % inicio
f=length(fail); % fim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot fasores da componente fundamental
figure
subplot(2,1,1), hold on, grid on
%plot(ref_seg(2,:),'o-')
loglog( amp_med(1,:), '*-')
loglog( amp_med(2,:), 'o-')
loglog( amp_med(3,:), 'd-')
%title('Amplitudes (pu)')
%legend('DFT','R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)')
ylabel('Amplitude (pu)'), xlabel('tempo (s)')
subplot(2,1,2), hold on, grid on

plot( (fase_med(1,:)), '*-')
plot( (fase_med(2,:)), 'o-')
plot( (fase_med(3,:)), 's-')
%title('Fase (rad)')
%legend('DFT','R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)')
ylabel('Fase (rad)'), xlabel('tempo (s)')


% get(gca,'fontname')  % shows you what you are using.
% set(gca,'fontname','times')  % Set it to times

% Plot ERRO de fasores (DFT e Prony)
figure  
subplot(3,1,1), hold on, grid on%
plot( tve(1,:), '*-')
plot( tve(2,:), 'o-')
plot( tve(3,:), 's-')
ylabel('TVE (%)'), xlabel('expoente')
%title('TVE (%)') 
%legend('DFT', 'R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)')
%title('TVE (%)'), legend('R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)') 
subplot(3,1,2), hold on, grid on
plot( amp_error(1,:), '*-')
plot( amp_error(2,:), 'o-')
plot( amp_error(3,:), 's-')
%title('Erro de amplitude (%)')
%legend('DFT', 'R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)')
%title('Erro de amplitude (%)'), legend('R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)') 
ylabel('Erro de amplitude (%)'), xlabel('expoente')
subplot(3,1,3), hold on, grid on
loglog( phase_error(1,:), '*-')
loglog( phase_error(2,:), 'o-')
loglog( phase_error(3,:), 's-')
%title('Erro de fase (º)')
legend('DFT', 'Sparse (OMP1)', 'Sparse (OMP2)')
%title('Erro de fase (º)'), legend('R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)') 
ylabel('Erro de fase (º)'), xlabel('expoente')



%x = logspace(1,20,20);
x = 1:1:20;
figure
subplot(2,1,1, 'XScale', 'log', 'YScale', 'log')
loglog(abs(phase_error(2,:)),  '*-')
title('Correção: 10^{-p}')
xlabel('expoente p'), ylabel('Erro de amplitude (%)')
%title(['Erro de amplitude: N=' num2str(N)])
grid on

subplot(2,1,2, 'XScale', 'log', 'YScale', 'log')
loglog(abs(amp_error(2,:)/100),  '*-')
%title('Erro de fase (º)')
xlabel('expoente p'), ylabel('erro de fase (º)')
grid on
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:2
   max_tve(i) = max(tve(i,:));
   max_amp_error(i) = max(amp_error(i,:)); 
   max_phi_error(i) = max(phase_error(i,:)); 
end
%