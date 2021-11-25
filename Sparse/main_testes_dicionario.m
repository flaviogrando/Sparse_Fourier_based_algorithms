%%% TESTES ESTIMADOR DFT-PRONY

clc
clear all
close all


set(0, 'defaultAxesFontSize',12);
set(0, 'defaultAxesFontName','times');

wb = waitbar(0,'Progresso');

ff=-5.12:0.01:5.11; % cria vetor para deslocar o ponto da grade

for w=1:length(ff)

waitbar(w/length(ff));
% ---------------------------------------------------------------------
% INICIALIZAÇÕES

N = 256;       % número de amostras
m = 1;       % passo de janelamento
Fs = 12800;    % Taxa de amostragem
f0 = 50;      % freq. fundamental (teórica)
noise = 00;    % nível de ruído em dB ( 0 = sem ruído)
T = 0.02+1/Fs;      % tempo total (em segundos)
param = ff(w);    % amort(w);    % freq. de modulação, por exemplo
type = 4;     % seleciona tipo de teste (ver gerador_sinais.m)

freqs = [50];
% freqs = [49.5 99 148.51 198.01];        % Vetor de frequências

%phases = zeros(3,2);
% phases = [3.05 4.85 4.75 3.11;           % vetor de fases
%          -2*pi/3 0 0 0;
%          +2*pi/3 0 0 0];    
phases(1,1) = deg2rad(0);
phases(2,1) = deg2rad(-120);
phases(3,1) = deg2rad(+120);

phases(1,2) = deg2rad(0);
phases(2,2) = deg2rad(-120);
phases(3,2) = deg2rad(+120);

amps = [1 0.1; 1 0.1; 1 0.1;];
%  amps = [6 2.21 0.812 0.299;           % vetor de amplitudes
%          1 0.3 0 0;
%          1 0.3 0 0]; 

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
eixo_x = 0:granN:granN*(N-1);
% 
% % Plot espectro
% figure
% hold on, grid on
% stem(eixo_x,2*abs(Sa(1,:))/N)
% % stem(abs(Sa))
% % stem(abs(Sa))
% title('Amplitude (pu) - DFT')
% ylabel('pu'), xlabel('Frequência (Hz)')


Sa = 2*(Sa)/N;

% Sa(1,:) = Sa(1,:)/sum(Sa(1,:));



Dic(:,w) = Sa(1,:).';


end
close(wb)
% 
% figure, grid on, stem(Dic(1,:))
% figure, grid on, stem(Dic(2,:))
% figure, grid on, stem(Dic(3,:))
% figure, grid on, stem(Dic(255,:))
% 
% figure, grid on, stem(imag(Dic(1,:)))
% figure, grid on, stem(imag(Dic(2,:)))
% figure, grid on, stem(imag(Dic(3,:)))
% figure, grid on, stem(imag(Dic(255,:)))
figure, grid on, stem(abs(Dic(1,:)))
figure, grid on, stem(abs(Dic(2,:)))
figure, grid on, stem(abs(Dic(3,:)))
figure, grid on, stem(abs(Dic(255,:)))

% figure, grid on, stem(imag(Dic(,:)))
% figure, grid on, stem(imag(Dic(2,:)))
% figure, grid on, stem(imag(Dic(3,:)))
% figure, grid on, stem(imag(Dic(255,:)))


figure, hold on, grid on
stem(Dic(:,3));
stem(Dic(:,1022));
% subplot(2,1,1), hold on, grid on, stem(Dic(:,256));
% subplot(2,1,2), hold on, grid on, stem(Dic(:,768));

figure
mesh(abs(Dic(3:255,:)))

% save('Dic_A.mat', 'Dic');
