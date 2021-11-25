%%%% FUN��O DICION�RIO


function [Dic ] = dic_time_adapt(N, H, Fs, f0, freq)

Dic = zeros(N,H);
norm1_col = zeros(H);
GranH = (Fs/(H));

Delta_f = freq - f0; % desvio de frequ�ncia
if Delta_f > GranH
    % Fra��o do grid - extrai parte fracion�ria
    delta = Delta_f/GranH - fix(Delta_f/GranH);
elseif Delta_f < -GranH
    delta = Delta_f/GranH - fix(Delta_f/GranH);
else
    delta = Delta_f/GranH;
end

%     % Menor deslocamento possivel
%     if delta(j) >= 0.5
%         delta(j) = -(1 - delta(j));
%     elseif delta(j) <-0.5
%         delta(j) = 1 + delta(j);
%     end

resol = 1/P;
k = 0:resol:N-resol;

k = k + delta;

i = 1;
for n=0:N-1;
    Dic(i,:) = exp(-1j*2*pi*k*n/N);
    %Dfm(i,:) = cos(2*pi*k*n/N)-1j*sin(2*pi*k*n/N);
    i=i+1;
end
% Normaliza��o
for i=1:length(k)
    norm1_col(i) = sum(abs(Dic(:,i)));
    Dic(:,i) = Dic(:,i)./(norm1_col(i));
end

end