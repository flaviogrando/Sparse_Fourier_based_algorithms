%%%% OMP
% implementação propria, ajustado para dicionário freq domain


function [amps, a, sup, res] = OMP_freq_freqA(k, signal, Dic, gran, freq)

[m,n] = size(Dic); % Get dictionary dimension

% 1) initializations
error_goal = 1e-3;
amps = zeros(1,n);
%a = zeros(1,m);
% sup = zeros(1,m);
D = [];

r = signal;   % residue
% norm = sum(abs(r));
% r = r./norm;

i=0;
j=0;



%while (error > error_goal) || (i<20),
while (j<k)
    
    i = i+1;
    j = j+1; % iteração do resíduo


     if i==100 % PRIMEIRA ITERAÇÃO - RECUPERA FUNDAMENTAL
           
        % 2) obtain support vector (determine indexes)    
        sup(1) = 10 + 1;        
        sup(2) = n - 10+1;
        
        % 3) construct matrix D
        D = [Dic(:,sup(1)) Dic(:,sup(2))];
        
        % 4) amplitude adjustment (least squares)
        a = pinv(D)*signal;
        
        % atualiza o vetor de amplitudes (não itera mais)
        amps(sup(1)) = a(1);
        amps(sup(2)) = a(2);
        
        i = i+1; % atualiza uma iteração do suporte
        
    else % DEMAIS ITERAÇÕES - RECUPERA HARMÔNICOS
        
      
        % 2) obtain support vector (search indexes)
        temp = Dic(1:m/2,1:n/2)'*r(1:m/2,j);
        [~, index] = max(abs(temp(1:n/2)));
        sup(i) = index(1);
        sup(i+1) = n - index(1) + 2; %( +1 pra compensar o index)
        
        if sup(i+1)>n
            sup(i+1) = n;
        end
        
        % 3) construct matrix D
        D = [Dic(:,sup(i)) Dic(:,sup(i+1))];
        
        % 4) amplitude adjustment (least squares)
        a = pinv(D)*r(:,j);

        amps(sup(i)) = a(1);
        amps(sup(i+1)) = a(2);
        
        i = i+1; % atualiza uma iteração

    end
    % 5) new residual calculation
    ajuste(:,j) = D*a;
    %r(:,i+1) = signal - ajuste(:,i);
    r(:,j+1) = r(:,j) - ajuste(:,j);

end

    
for i=1:length(sup)
    %amps(sup(i)) = a(i);
    a(i) = amps(sup(i)); 
end

res = r;
% sup = sup;

end