%%%% função estimação esparsa com OMP
% Entrada: Spectro (DFT) - fases a, b e c
% Dic -  dicionário
% K - esparsidade
% P -  fator de interpolação (para o dicionário)

%function [Ssa, Ssb, Ssc, Asa, Asb, Asc, phisa, phisb, phisc, indx_a] = estimador_freq_omp(Sa, Sb, Sc, K, P, omp)
function [Ssa, Asa, phisa, indx_a] = estimador_freq_omp(Sa, K, P, omp, Fs, freq, Dic)

[num_win, win_size] = size(Sa);
N = win_size;
H = N*P;
granH = Fs/H;

% Cria variáveis
Ssa = zeros(num_win,H);
%Ssb = zeros(num_win,H);
%Ssc = zeros(num_win,H);
% sup_a = zeros(num_win,K);
% sup_b = zeros(num_win,K);
% sup_c = zeros(num_win,K);


mean_a = zeros(num_win,1);  % 1 = só a fundamental
%mean_b = zeros(num_win,1);
%mean_c = zeros(num_win,1);


%phisb = zeros(1,num_win);
%phisc = zeros(1,num_win);
    
% Cria dicionário DIRICHLET KERNEL
%[Dic] = dic_freq(N, H);
%load('E:\block_dic32u.mat') % carrega dicionário
%load('E:\block_dic256u.mat') % carrega dicionário

% Aplica OMP em cada espectro i de cada fase (a, b ou c)
switch omp  % Seleciona  OMP
    case 1% (Phi,V,m,alpha,errorGoal)
        for i=1:num_win
            [Ssa(i,:), coef(i,:), indx_a(i,:), ra] = OMP_freq_ord(K, Sa(i,:).', Dic);
            %[Ssa(i,:), indx_a(i,:),~,~, coef(i,:)] = OMP_1(Dic, Sa(i,:).',K,0,1e-10);
            %[Ssb(i,:), indx_b(i,:),~,~, coef_b(i,:)] = OMP_1(Dic, Sb(i,:).',K,0,1e-10);
            %[Ssc(i,:), indx_c(i,:),~,~, coef_c(i,:)] = OMP_1(Dic, Sc(i,:).',K,0,1e-10);
        end    
    case 2  
        for i=1:num_win
            [Ssa(i,:), coef(i,:), indx_a(i,:), ra] = OMP_freq_ord2(K, Sa(i,:).', Dic);
            %[Ssa(i,:), indx_a(i,:)] = OMP_2(K, Sa(i,:), Dic);
            %[Ssb(i,:), indx_b(i,:)] = OMP_2(K, Sb(i,:), Dic);
            %[Ssc(i,:), indx_c(i,:)] = OMP_2(K, Sc(i,:), Dic);
        end
    case 3
        %tic;
        for i=1:num_win
            [Ssa(i,:), coef(i,:), indx_a(i,:), ra] = OMP_freq(K, Sa(i,:).', Dic);
            %[Ssb(i,:), indx_b(i,:), rb] = OMP_freq2(K, Sb(i,:).', Dic);
            %[Ssc(i,:), indx_c(i,:), rc] = OMP_freq2(K, Sc(i,:).', Dic);           
        end
        %toc;
    case 4 
        for i=1:num_win
             [Ssa(i,:), coef(i,:), indx_a(i,:)] = OMP_freq_freq(K, Sa(i,:).', Dic, granH, freq(i));
            %[Ssb(i,:), indx_b(i,:), rb] = OMP_freq2(K, Sb(i,:).', Dic);
            %[Ssc(i,:), indx_c(i,:), rc] = OMP_freq2(K, Sc(i,:).', Dic);           
        end 
end

sh = length( coef(1,:));
Asa = zeros(num_win,sh);
phisa = zeros(num_win,sh);

for i=1:num_win
%     % Extrai a média das componentes de frequência espelhadas
%     mean_a(i,:) = (Ssa(i,sup_a(i,1)) + Ssa(i,sup_a(i,2)) )/2;
%     mean_b(i,:) = (Ssb(i,sup_b(i,1)) + Ssb(i,sup_b(i,2)) )/2;
%     mean_c(i,:) = (Ssc(i,sup_c(i,1)) + Ssc(i,sup_c(i,2)) )/2;
        
    % pega o valor da fundamental
    %indice = min(indx_a(i,:)
    %mean_a(i,:) = (Ssa(i,min(indx_a(i,:))) );
    %mean_b(i,:) = (Ssb(i,min(indx_b(i,:))) );
    %mean_c(i,:) = (Ssc(i,min(indx_c(i,:))) );
     
    % Calcula amplitude
    Asa(i,:) = abs(coef(i,:));
    %Asb(:,i) = abs(mean_b(i,:));
    %Asc(:,i) = abs(mean_c(i,:));
    
    % Calcula fase
    phisa(i,:) = angle(coef(i,:));
    %phisb(:,i) = angle(mean_b(i,:));
    %phisc(:,i) = angle(mean_c(i,:));
end



end