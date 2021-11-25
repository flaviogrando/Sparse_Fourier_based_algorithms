%%%% fun��o estima��o esparsa com OMP
% Entrada: Spectro (DFT) - fases a, b e c
% Dic -  dicion�rio
% K - esparsidade
% P -  fator de interpola��o (para o dicion�rio)

%function [Ssa, Ssb, Ssc, Asa, Asb, Asc, phisa, phisb, phisc, indx_a] = estimador_freq_omp(Sa, Sb, Sc, K, P, omp)
function [Ssa, Asa, phisa, indx_a] = estimador_time_omp(Sa, K, P, omp)

[num_win, win_size] = size(Sa);
N = win_size;
H = N*P;

% Cria vari�veis
Ssa = zeros(num_win,H);
%Ssb = zeros(num_win,H);
%Ssc = zeros(num_win,H);
% sup_a = zeros(num_win,K);
% sup_b = zeros(num_win,K);
% sup_c = zeros(num_win,K);


mean_a = zeros(num_win,1);  % 1 = s� a fundamental
%mean_b = zeros(num_win,1);
%mean_c = zeros(num_win,1);

Asa = zeros(1,num_win);
%Asb = zeros(1,num_win); 
%Asc = zeros(1,num_win); 

phisa = zeros(1,num_win);
%phisb = zeros(1,num_win);
%phisc = zeros(1,num_win);
    
% Cria dicion�rio DIRICHLET KERNEL
[Dic] = dic_time(N, H, P);


% Aplica OMP em cada espectro i de cada fase (a, b ou c)
switch omp  % Seleciona  OMP
    case 1% (Phi,V,m,alpha,errorGoal)
        for i=1:num_win
            [Ssa(i,:), indx_a(i,:),~,~, coef_a(i,:)] = OMP_1(Dic, Sa(i,:).',K,0,1e-10);
            %[Ssb(i,:), indx_b(i,:),~,~, coef_b(i,:)] = OMP_1(Dic, Sb(i,:).',K,0,1e-10);
            %[Ssc(i,:), indx_c(i,:),~,~, coef_c(i,:)] = OMP_1(Dic, Sc(i,:).',K,0,1e-10);
        end    
    case 2  
        for i=1:num_win
            [Ssa(i,:), indx_a(i,:)] = OMP_2(K, Sa(i,:), Dic);
            %[Ssb(i,:), indx_b(i,:)] = OMP_2(K, Sb(i,:), Dic);
            %[Ssc(i,:), indx_c(i,:)] = OMP_2(K, Sc(i,:), Dic);
        end
    case 3 
        for i=1:num_win
            [Ssa(i,:), indx_a(i,:), ra] = OMP_freq(K, Sa(i,:).', Dic);
            %[Ssb(i,:), indx_b(i,:), rb] = OMP_freq2(K, Sb(i,:).', Dic);
            %[Ssc(i,:), indx_c(i,:), rc] = OMP_freq2(K, Sc(i,:).', Dic);           
        end       
end



for i=1:num_win
%     % Extrai a m�dia das componentes de frequ�ncia espelhadas
%     mean_a(i,:) = (Ssa(i,sup_a(i,1)) + Ssa(i,sup_a(i,2)) )/2;
%     mean_b(i,:) = (Ssb(i,sup_b(i,1)) + Ssb(i,sup_b(i,2)) )/2;
%     mean_c(i,:) = (Ssc(i,sup_c(i,1)) + Ssc(i,sup_c(i,2)) )/2;
        
    % pega o valor da fundamental
    mean_a(i,:) = (Ssa(i,min(indx_a(i,:))) );
    %mean_b(i,:) = (Ssb(i,min(indx_b(i,:))) );
    %mean_c(i,:) = (Ssc(i,min(indx_c(i,:))) );
     
    % Calcula amplitude
    Asa(:,i) = abs(mean_a(i,:));
    %Asb(:,i) = abs(mean_b(i,:));
    %Asc(:,i) = abs(mean_c(i,:));
    
    % Calcula fase
    phisa(:,i) = angle(mean_a(i,:));
    %phisb(:,i) = angle(mean_b(i,:));
    %phisc(:,i) = angle(mean_c(i,:));
end



end