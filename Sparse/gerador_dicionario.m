%%% GERADOR DE DICIONÁRIO
clear all
close all

wb = waitbar(0,'Progresso');


N = 256;
P = 10;
H = N*P;
Fs = 12800;
GranH = (Fs/(H));
P2 = 50;
GranH2 = Fs/(P2*H);

delta = -0.5:(1/P2):0.5;

Dic = zeros(N,H);

ajuste = (1e-12);
block = [];

%%%%%%%%%% BLOCO DE DICIONÁRIOS PELO AJUSTE DELTA (RUIM)
% for i=1:P2+1;
%     
%     waitbar(i/P2);
%     
%     for l=0:H-1      % colunas
%         for k=0:N-1  % linhas             
%             arg = (k/N - (l+ajuste+delta(i))/H); % argumento do kernel 
%             Dic(k+1,l+1) = (sin(pi*N*(arg))/(N*sin(pi*(arg))))*exp(-1i*pi*(N-1)*(arg)); % eq. (3)-1
%             %Dic(k+1,l+1) = (sin(pi*N*(arg))/(N*sin(pi*(arg))));
%             %Dic(k+1,l+1) = sinc(arg);
%     %         if isnan(Dic(k+1,l+1))
%     %             Dic(k+1,l+1) = 1;
%     %         end     
%         end
%     end
%     block(:,:,i) = Dic;
% %     field = ['dic',num2str(i)];
% %     block = struct(field,Dic);
% end
% close(wb);
% save('E:\block_dic32m.mat', 'block');
 
 
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
        %Dic(k+1,l+1) = (sin(pi*N*(arg))/(N*sin(pi*(arg))));
        %Dic(k+1,l+1) = sinc(arg);
        grade_arg(l+1) = arg;
        
        if isnan(Dic(k+1,l+1))
            Dic(k+1,l+1) = (-1)^(k*(N-1));
        end     
    end
end
close(wb);
save('D:\dicionarios_cs_dft\block_dic256u.mat', 'Dic');
figure
mesh(abs(Dic))



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

% load('E:\block_dic256u.mat')
% 
 %%%%%%%%%%%%%%%%%% BLOCO DE DICIONÁRIOS 2
N = 256;
P = 10;
P2 = 100;
H = N*P;
H2 = N*P*P2;

granH = (Fs/(H)); % Nova granularidade
gradeH = 0:granH:granH*((H)-1);
granH2 = (Fs/(H2)); % Nova granularidade
gradeH2 = 0:granH2:granH2*((H2)-1);

num_dic = 200;
for n=0:num_dic;
    
    for l=0:H/2-1;
        ln = l*P2 + n;
        
        gradeH2_map(n+1,l+1) = ln*granH2;
        
        new_Dic(:,l+1) = Dic(:,ln+1);
        new_Dic(:,H-l) = Dic(:,H2 - P2 + 1 - ln);
    end
    block(:,:,n+1) = new_Dic;
end

save('D:\dicionarios_cs_dft\block_dic256us.mat', 'block');

