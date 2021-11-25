%%% GERADOR DE DICIONÁRIO
clear all
close all

wb = waitbar(0,'Progresso');


N = 256;          % Número de amostras
P = 10;           % Fator de interpolação
H = N*P;          % Nova dimensão
Fs = 12800;       % Taxa de amostragem
GranH = (Fs/(H)); % Nova granularidade
P2 = 50;
GranH2 = Fs/(P2*H);

delta = -0.5:(1/P2):0.5;

Dic = zeros(N,H);

block = [];
 
 
%%%%%%%%%%%%%%%%% DICIONÁRIO ÚNICO
N = 256;
P = 1000;
H = N*P;

Dic = zeros(N,H);
ajuste = (1e-9);

for k=0:N-1      % colunas
    
    waitbar(k/N);
    
    for l=0:H-1  % linhas
        
        arg = (k/N - (l+ajuste)/H); % argumento do kernel 
        Dic(k+1,l+1) = (sin(pi*N*(arg))/(N*sin(pi*(arg))))*exp(-1i*pi*(N-1)*(arg)); % eq. (3)-1
        grade_arg(l+1) = arg;
        
        % Inserindo 1 no lugar de NaN
        if isnan(Dic(k+1,l+1))
            Dic(k+1,l+1) = (-1)^(k*(N-1));
        end     
    end
end
close(wb);
save('C:\dicionarios_cs_dft\block_dic256u.mat', 'Dic');
% figure
% mesh(abs(Dic))

% 
% figure
% mesh(abs(block(:,:,1)))
% figure
% mesh(abs(block(:,:,50)))
% 
% figure, grid on, hold on
% stem(abs(block(:,10,1)))   % 45
% stem(abs(block(:,10,51)))  % 50
% stem(abs(block(:,10,101))) % 55
% %subplot(2,2,1), grid on, hold on
% 
% figure, grid on, hold on
% stem(abs(block(:,320-10+1,1)))   % 45
% stem(abs(block(:,320-10+1,51)))  % 50
% stem(abs(block(:,320-10+1,101))) % 55




 %%%%%%%%%%%%%%%%%% BLOCO DE DICIONÁRIOS - Dicionários Múltiplos
    
N = 256;             % Taxa de amostragem  
P = 10;              % Fator de interpolação (intra-dicionário)
P2 = 100;            % Fator de interpolação (inter-dicionário)
H = N*P;             % Granularidade intra-dicionário
H2 = N*P*P2;         % Granularidade inter-dicionário

num_dic = 200;       % Número de dicionários

granH = (Fs/(H)); % Nova granularidade
gradeH = 0:granH:granH*((H)-1);
granH2 = (Fs/(H2)); % Nova granularidade
gradeH2 = 0:granH2:granH2*((H2)-1);


for n=0:num_dic
    
    for l=0:H/2-1
        ln = l*P2 + n;
        
        gradeH2_map(n+1,l+1) = ln*granH2;
        
        new_Dic(:,l+1) = Dic(:,ln+1);
        new_Dic(:,H-l) = Dic(:,H2 - P2 + 1 - ln);
    end
    block(:,:,n+1) = new_Dic;
end

save('C:\dicionarios_cs_dft\block_dic256us.mat', 'block');

