%%%% FUNÇÃO DICIONÁRIO


function [Dic] = dic_freq_seletor(N, H, Fs, f0, freq, block)

%load('E:\block_dic.mat') % carrega dicionários

% Dic = zeros(N,H);
GranH = (Fs/(H));
P2 = 100;
granH2 = Fs/(P2*H);

% freq_mtx = ones(P2,3)*freq;
% 
% Delta_f = freq_mtx  - gradeH2(:,10:12); 
% 
% 
% [colunas,linhas] = min(abs(Delta_f));
% 
% [~, min_col_indx] = min(colunas);
% 
% atomo = 9 + min_col_indx;
% 
% selection = linhas(min_col_indx)
% 
% freq_selected = gradeH2(selection,atomo);

%[valor, indice] = min(Delta_f(:,indx(1)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vetor_freq = 45:0.05:55;

delta_fn = vetor_freq - freq;

[~,selection] = min(abs(delta_fn));

% vetor_freq = -GranH/2:(GranH/P2):GranH/2;
% 
% Delta_f = freq - f0; % desvio de frequência
% 
% delta_fn = vetor_freq - Delta_f;%/GranH;
% 
% [~,selection] = min(abs(delta_fn))
% delta_fn = round(Delta_f/GranH2);
% 
% if delta_fn >=0
%     selection = delta_fn+1;
% elseif delta_fn < 0
%     selection = 10+delta_fn+1;
% end

Dic = block(:,:,selection);


end



