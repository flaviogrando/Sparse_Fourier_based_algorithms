%%%% OMP
% implementação propria, ajustado para dicionário freq domain


function [amps, a, sup, res] = OMP_freq_orig(k, signal, Dic)

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
while (i<k)
    
    i = i+1;
    %b 2) obtain support vector (search indexes)
    temp = Dic'*r(:,i);
    [~, index] = max(abs(temp));
    sup(i) = index(1);
    
 
    % 3) construct matrix D
    d = Dic(:,sup(i));
    D = [D d];
    

    % 4) amplitude adjustment (least squares)
    a = pinv(D)*signal;

    % 5) residue 
    ajuste(:,i) = D*a;
    r(:,i+1) = signal - ajuste(:,i);

end

     
for i=1:k;
    
    amps(sup(i)) = a(i);
        
end

res = r;
% sup = sup;

end