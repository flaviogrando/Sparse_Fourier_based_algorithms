%%%% FUNÇÃO DICIONÁRIO


function [Dic] = dic_freq_adapt(N, H, Fs, f0, freq)

Dic = zeros(N,H);

P = 10;
P2 = 50;
H = N*P;
H2 = N*P*P2;

GranH = (Fs/(H));
GranH2 = (Fs/(H2));
Delta_f = freq - f0; % desvio de frequência
if Delta_f > GranH
    % Fração do grid - extrai parte fracionária
    delta = Delta_f/GranH - fix(Delta_f/GranH);
elseif Delta_f < -GranH
    delta = Delta_f/GranH - fix(Delta_f/GranH);
else
    delta = Delta_f/GranH;
end

% % Menor deslocamento possivel
% if delta >= 0.5
%     delta = -(1 - delta);
% elseif delta <-0.5
%     delta = 1 + delta;
% end

ajuste = 1e-12; %(1e-12);
for k=0:N-1      % colunas
    for l=0:H/2-1  % linhas             
        arg = (k)/(N) - ((l + ajuste + delta)/(H) ); % argumento do kernel
        Dic(k+1,l+1) = (sin(pi*N*(arg))/(N*sin(pi*(arg))))*exp(-1i*pi*(N-1)*(arg)); % eq. (3)-1                        
    end
    for l=H/2:H-1  % linhas             
        arg = (k)/(N) - ((l + ajuste - delta)/(H) ); % argumento do kernel
        Dic(k+1,l+1) = (sin(pi*N*(arg))/(N*sin(pi*(arg))))*exp(-1i*pi*(N-1)*(arg)); % eq. (3)-1                        
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ajuste = 1e-9;
% 
% vetor_freq = 45:0.1:55;
% 
% delta_f = vetor_freq - freq;
% 
% [~,n] = min(abs(delta_f));
% 
% delta = -delta_f(n)*(1/GranH2);
% 
% for k=0:N-1      % linhas
%     
%     for l=0:H/2-1  % linhas 
%         
%         ln = l*P2 + n - 51;   % índice refinado com passos P2;
%         
%         %gradeH2_map(l+1) = ln*GranH2;
%         
%         arg = (k)/(N) - ((ln + ajuste + delta )/(H2) ); % argumento do kernel
%         Dic(k+1,l+1) = (sin(pi*N*(arg))/(N*sin(pi*(arg))))*exp(-1i*pi*(N-1)*(arg));
%         
%        % grade_arg(l+1) = arg;
%         %ln = H2 - P2 + 1 - (l*P2 + n);
%         
%         arg = (N-k)/(N) - (H2 + P2*P - ln + ajuste - P2 - delta)/(H2) ; % argumento do kernel
% 
%         %arg = 1 - arg;
%         Dic(N-k,H-l) = (sin(pi*N*(arg))/(N*sin(pi*(arg))))*exp(-1i*pi*(N-1)*(arg));
%         
%         %grade_arg(H - l) = arg;
%     end
% end


% for k=0:N-1      % linhas
%     
%     for l=0:H/2-1  % linhas 
% 
%         ln = l*P2 + n - 50 - 1;   % índice refinado com passos P2;
%         
%         gradeH2_map(l+1) = ln*GranH2;
%         
%         arg = (k)/(N) - ((ln + ajuste + delta )/(H2) ); % argumento do kernel
%         Dic(k+1,l+1) = (sin(pi*N*(arg))/(N*sin(pi*(arg))))*exp(-1i*pi*(N-1)*(arg));
%           
%     end
%     
% %     for l=0:H/2-1
% %         Dic(k+1,l+1)
% %     end
%     
% end


% figure
% mesh(abs(Dic))
% 
% figure
% mesh(angle(Dic))
% 
% display(n)

% num_dic = 100;
% for n=0:num_dic;
%     
%     for l=0:H/2-1;
%         ln = l*P2 + n;
%         
%         gradeH2_map(n+1,l+1) = ln*granH2;
%         
%         new_Dic(:,l+1) = Dic(:,ln+1);
%          Dic(:,l+1) = Dic(:,ln+1);
%         new_Dic(:,H-l) = Dic(:,H2 - P2 + 1 - ln);
%     end
%     block(:,:,n+1) = new_Dic;
% end