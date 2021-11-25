%%%% OMP
% implementação propria, ajustado para dicionário time domain


function [amps, sup, res] = OMP_freq(k, signal, Dic);

[m,n] = size(Dic); % Get dictionary dimension

% 1) initializations
error_goal = 1e-3;
amps = zeros(1,n);
%a = zeros(1,m);
% sup = zeros(1,m);
D = [];

r = signal;   % residue
error = sum(abs(r));

i=0;

%while (error > error_goal) || (i<20),
while (i<k)
    
    i = i+1;
    % 2) obtain support vector (search indexes)
    temp = Dic.'*r(:,i);
    [~, index] = max(abs(temp));
%     if i==1
%         [~, index] = max(abs(temp(n/2:n))); % seleciona somente metade do espectro
%     else
%         [~, index] = max(abs(temp(1:n/2))); % seleciona somente metade do espectro
%     end
    sup(i) = index(1);

%     figure
%     subplot(3,1,1), hold on, grid on
%     stem(r(:,i))
%     subplot(3,1,2), hold on, grid on
%     stem(Dic(:,10))
%     subplot(3,1,3), hold on, grid on
%     stem(Dic(:,12))
%     
%     figure
%     stem((temp))
    
    % 3) construct matrix D
    if sup(i)>0
        d = Dic(:,sup(i));
        D = [D d];
    end

    % 4) amplitude adjustment (least squares)
    a = pinv(D)*signal;
    %a = ((D.'*D)^(-1))*D.'*signal;
    %a = D.'*signal;

    % 5) new residual calculation
    r(:,i+1) = signal - D*a;
    error = sum(abs(r(i+1)));
    
%     figure
%     subplot(2,1,1), grid on
%     stem(r(:,1))
%     subplot(2,1,2), grid on
%     stem(r(:,2))
end

%display(i)

for i=1:length(a)
    amps(sup(i)) = a(i);
    
end

res = r;
% sup = sup;

end