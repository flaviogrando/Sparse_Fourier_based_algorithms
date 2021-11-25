%%%% OMP
% implementação propria, ajustado para dicionário freq domain


function [amps, a, sup, res] = OMP_freq(k, signal, Dic)

[m,n] = size(Dic); % Get dictionary dimension

% 1) initializations
error_goal = 1e-3;
amps = zeros(1,n);
%a = zeros(k*2,1);
sup = zeros(1,k*2);
%ajuste = zeros(m,k);
D = [];

r = signal;   % residue
% norm = sum(abs(r));
% r = r./norm;
% r = r.^2;
% Dic = Dic.^2;

i=0;
j=0;
%erro = sum(abs(r));

%while erro > 0.1;
while (j<k)
    
    i = i+1;
    j = j+1;

    % 2) obtain support vector (search indexes)
    temp = abs(Dic(1:m/2,1:n/2))'*abs(r(1:m/2,j));
    [~, index] = max(abs(temp(1:n/2)));
    sup(i) = index(1);
    sup(i+1) = n - index(1) + 2; %( +1 pra compensar o index)

    % 3) construct matrix D
    %D = [Dic(:,sup(i)), Dic(:,sup(i+1))];
    D = [D Dic(:,sup(i)) Dic(:,sup(i+1))];

    % 4) amplitude adjustment (least squares)
    %ah = pinv(D)*r(:,j);
    ah = pinv(D)*signal;

%         amps(sup(i)) = ah(1);
%         amps(sup(i+1)) = ah(2);
    amps(sup(i)) = ah(i);
    amps(sup(i+1)) = ah(i+1);

    % atualiza 'a' com as harmônicas, sem interferir na fundamental
    a(i,:) = ah(i);
    a(i+1,:) = ah(i+1);
% atualização refinada a excluindo a funtamental (utlizado em caso de
% harmonicos
    for l=3:length(ah)
        a(l) = ah(l);
    end 

   i = i+1; % atualiza uma iteração (em caso de um par de freq por vez)
       


    % 5) new residual calculation
    ajuste(:,j) = D*a;
    r(:,j+1) = signal - ajuste(:,j);


end

     
for i=1:k;
    
    amps(sup(i)) = a(i);
        
end

res = r;
% sup = sup;

end