%%%% FUNÇÃO DICIONÁRIO


function [Dic ] = dic_freq(N, H)

Dic = zeros(N,H);

ajuste = (1e-9);

for l=0:H-1      % colunas
    for k=0:N-1  % linhas             
        arg = (k/N - (l+ajuste)/H); % argumento do kernel 
        %arg = arg*(2*pi);
        Dic(k+1,l+1) = (sin(pi*N*(arg))/(N*sin(pi*(arg))))*exp(-1i*pi*(N-1)*(arg)); % eq. (3)-1
        %Dic(k+1,l+1) = Dic(k+1,l+1)*1/N;
        %Dic(k+1,l+1) = (sin(pi*N*(arg))/(N*sin(pi*(arg))));
        %Dic(k+1,l+1) = diric(arg,N)*exp(-1i*pi*(N-1)*(arg));
        %Dic(k+1,l+1) = sinc(arg);
%         if isnan(Dic(k+1,l+1))
%             Dic(k+1,l+1) = 1;
%         end     
    end


%load('Dic_A.mat')
end

 %save('Dic_freq.mat', 'Dic');