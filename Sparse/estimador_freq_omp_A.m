%%%% função estimação esparsa com OMP
% Entrada: Spectro (DFT) - fases a, b e c
% Dic -  dicionário
% K - esparsidade
% P -  fator de interpolação (para o dicionário)

function [Ssa, Asa, phisa, indx_a] = estimador_freq_omp_A(Sa, K, P, omp, Fs, f0, freq)

[num_win, win_size] = size(Sa);
N = win_size;
H = N*P;
granH = (Fs/(P*N));


% Cria variáveis
Ssa = zeros(num_win,H);
%Ssb = zeros(num_win,H);
%Ssc = zeros(num_win,H);
% sup_a = zeros(num_win,K);
% sup_b = zeros(num_win,K);
% sup_c = zeros(num_win,K);
%Dic = zeros(N,H);


    
% w = flattopwin(H);

%for j=1:length(freq)
    
% Cria dicionário DIRICHLET KERNEL adaptado a frequencia
%[Dic] = dic_freq_adapt(N, H, Fs, f0, freq(j));    


    % Aplica OMP em cada espectro i de cada fase (a, b ou c)
    switch omp  % Seleciona  OMP
        case 1% (Phi,V,m,alpha,errorGoal)
            for i=1:num_win
                [Ssa(i,:), indx_a(i,:),~,~, coef_a(i,:)] = OMP_1(Dic, Sa(i,:).',K,0,1e-4);
                %[Ssb(i,:), indx_b(i,:),~,~, coef_b(i,:)] = OMP_1(Dic, Sb(i,:).',K,0,1e-4);
                %[Ssc(i,:), indx_c(i,:),~,~, coef_c(i,:)] = OMP_1(Dic, Sc(i,:).',K,0,1e-4);
            end    
        case 2  
            for i=1:num_win
                [Dic] = dic_freq_adapt(N, H, Fs, f0, freq(i));           % <<<< CONSTRÓI DICIONÁRIO
                [Ssa(i,:), coef(i,:), indx_a(i,:)] = OMP_freq(K, Sa(i,:).', Dic);
                %[Ssb(i,:), indx_b(i,:)] = OMP_2(K, Sb(i,:), Dic);
                %[Ssc(i,:), indx_c(i,:)] = OMP_2(K, Sc(i,:), Dic);
            end
        case 3
            %tic;
            for i=1:num_win     % OMP tradicional
                [Dic] = dic_freq_adapt(N, H, Fs, f0, freq(i));           % <<<< CONSTRÓI DICIONÁRIO
                [Ssa(i,:), coef(i,:), indx_a(i,:)] = OMP_freq(K, Sa(i,:).', Dic);
                %[Ssb(i,:), indx_b(i,:), rb] = OMP_freq(K, Sb(i,:).', Dic);
                %[Ssc(i,:), indx_c(i,:), rc] = OMP_freq(K, Sc(i,:).', Dic);               
            end
            %toc;
        case 4
            for i=1:num_win     % OMP com seleção de atomo forçada
                [Dic] = dic_freq_adapt(N, H, Fs, f0, freq(i));           % <<<< CONSTRÓI DICIONÁRIO
                [Ssa(i,:), coef(i,:), indx_a(i,:)] = OMP_freq_freqA(K, Sa(i,:).', Dic, granH, freq(i));
            end
        case 5
            for i=1:num_win
                [Ssa(i,:), indx_a(i,:)] = OMP_freq2(K, Sa(i,:).', Dic, GranH, P, freq(i));
            end
    end

%end
sh = length( coef(1,:));
Asa = zeros(num_win,sh);
phisa = zeros(num_win,sh);

for i=1:num_win
%     % Extrai a média das componentes de frequência espelhadas
%     mean_a(i,:) = (Ssa(i,sup_a(i,1)) + Ssa(i,sup_a(i,2)) )/2;
%     mean_b(i,:) = (Ssb(i,sup_b(i,1)) + Ssb(i,sup_b(i,2)) )/2;
%     mean_c(i,:) = (Ssc(i,sup_c(i,1)) + Ssc(i,sup_c(i,2)) )/2;
    
    % atualiza espectro Ssa com valor 'coef_a' no índice 'indx_a'
%     for j=1:n;
%         Ssa(i,indx_a(i,j)) = coef_a(i,j); 
%     end
    
    % Calcula amplitude
    Asa(i,:) = abs(coef(i,:));
    %Asb(:,i) = abs(mean_b(i,:));
    %Asc(:,i) = abs(mean_c(i,:));
    
    % Calcula fase
    phisa(i,:) = angle(coef(i,:));
    %Asb(:,i) = abs(mean_b(i,:));
    %Asc(:,i) = abs(mean_c(i,:));
    

end



end