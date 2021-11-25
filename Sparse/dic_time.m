%%%% FUNÇÃO DICIONÁRIO


function [Dic ] = dic_time(N, H, P)

Dic = zeros(N,H);

resol = 1/P;
k = 0:resol:N-resol;

i = 1;
for n=0:N-1;
    Dic(i,:) = exp(-1j*2*pi*k*n/N);
    %Dfm(i,:) = cos(2*pi*k*n/N)-1j*sin(2*pi*k*n/N);
    i=i+1;
end
% Normalização
for i=1:length(k)
    norm1_col(i) = sum(abs(Dic(:,i)));
    Dic(:,i) = Dic(:,i)./(norm1_col(i));
end

end