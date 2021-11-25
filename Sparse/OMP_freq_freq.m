%%%% OMP
% implementação propria, ajustado para dicionário freq domain


function [amps, a, sup, res] = OMP_freq_freq(k, signal, Dic, gran, freq)

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

    
    if i==1 % PRIMEIRA ITERAÇÃO - RECUPERA FUNDAMENTAL
     
        sup(1) = round(freq/gran) + 1;                   
        sup(2) = n- round(freq/gran) + 1;

        % 3) construct matrix D
        D = [Dic(:,sup(1)) Dic(:,sup(2))];
        
        % 4) amplitude adjustment (least squares)
        a = pinv(D)*signal;
        
        % Atualiza o vetor de amplitudes
        amps(sup(1)) = a(1);
        amps(sup(2)) = a(2);
        
        i = i+1; % atualiza uma iteração
        
    else % DEMAIS ITERAÇÕES - RECUPERA HARMÔNICOS
        
%         isodd = mod(i,2);
%         
%         % 2) obtain support vector (search indexes)
%         if isodd==1 % 1/2 quadrante (freq positivas)
%             temp = Dic(1:m/2,1:n/2)'*r(1:m/2,j);
%             [~, index] = max(abs(temp(1:n/2)));
%             sup(i) = index(1);
%         else     % 3/4 quadrante (freq positivas)% 
% %             temp = Dic(1:end,n/2:end)'*r(1:end,i);
% %             [~, index] = max(abs(temp));
% %             sup(i) = index(1)+n/2;
%             sup(i) = n - index(1) + 2;
%         end
%         
%         % 3) construct matrix D
%         d = Dic(:,sup(i));
%         D = [D d];
% 
%         % 4) amplitude adjustment (least squares)
%         if isodd==1
%             ah = pinv(D(1:m/2,:))*signal(1:m/2);
% %            amps(sup(i)) = ah(i); % atualiza as amplitudes
%             a(i) = ah(i);
%         else
%             ah = pinv(D)*signal;
% %            amps(sup(i-1)) = ah(i-1); % atualiza as amplitudes
% %            amps(sup(i)) = ah(i); % atualiza as amplitudes
% %            a(i-1) = ah(i-1);
%             a(i) = ah(i);
%         end
%         for l=3:length(ah)
%             a(l) = ah(l);
%         end 

        
        % 2) obtain support vector (search indexes)
        temp = Dic(1:m/2,1:n/2)'*r(1:m/2,j);
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
        
%         % atualiza 'a' com as harmônicas, sem interferir na fundamental
%         a(i) = ah(i);
%         a(i+1) = ah(i+1);
        for l=1:length(ah)
            a(l) = ah(l);
        end 
       
       i = i+1; % atualiza uma iteração (em caso de um par de freq por vez)
       
    end

    % 5) new residual calculation
    ajuste(:,j) = D*a;
    r(:,j+1) = signal - ajuste(:,j);
    %r(:,j+1) = r(:,j) - ajuste(:,j);
       
end

% for i=3:length(sup)
%     %amps(sup(i)) = a(i);
%     a(i) = amps(sup(i)); 
% end
    


res = r;
% sup = sup;

end