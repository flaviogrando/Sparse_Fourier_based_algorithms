%%%% OMP
% implementação propria, ajustado para dicionário freq domain

% MODIFICAÇÃO: TESTE DE MSE NO RESÍDUO DEPOIS DO AJUSTE DE AMPLITUDE
% PARA DIFERENTES ÁTOMOS

function [amps, sup, r] = OMP_freq(k, signal, Dic);

[m,n] = size(Dic); % Get dictionary dimension

P = 10; % AVALIAR CRITÉRIO PARA DETERMINAÇÃO DA MARGEM

% 1) initializations
error_goal = 1e-3;
amps = zeros(1,n);
%a = zeros(1,m);
sup = zeros(1,k);       % vetor de suporte
margin = floor(P/2);    % determina margem de índices admissível
r1  = ones(m,margin);   % vetor de residuos
mse_r1 = ones(1,margin);% vetor do MSE para o resíduo

% Inicializações
r = signal;             % residue
D = [];

i=0;

%while (error > error_goal) || (i<20),
while (i<k)
    
    i = i+1;
    % 2) obtain support vector (search indexes)
    %temp = Dic.'*r(:,i);
    temp = Dic.'*r;
    [~, index] = max(abs(temp));
    sup(i) = index(1);
    

    
    % Itera -margin/2 ate margin/2 átomos do escolhido anteriormente (index);
    lim = floor(margin/2)+1;             % define limite para inicio da iteração
    %figure
    for j=1:margin;

        index1(j) = sup(i)-lim+j;
 
        r0(:,j) = r - Dic(:,index1(j));   % calcula residuo para cada átomo
        mse_r0(:,j) = mean(abs(r0(:,j)).^2);% calcula MSE de cada resíduo
        
        %subplot(5,1,j), hold on, grid on, stem(Dic(:,index1(j))), title(['indice: ', num2str(index1(j))])
    
        a1(:,j) = pinv(Dic(:,index1(j)))*signal;
        
        r1(:,j) = r - Dic(:,index1(j))*a1(:,j);   % calcula residuo para cada átomo
        
        mse_r1(:,j) =  mean(abs(r1(:,j)).^2);% calcula MSE de cada resíduo
        
    end
    
%     figure, grid on, hold on, stem(abs(temp))
%     figure, stem(mse_r1), title(['K = ', num2str(i)])

    
    % Seleciona novo atomo com base no menor MSE
    [~, index] = min(mse_r1);
    sup(i) = index(1)+ index1(1)-1;

    
    % 3) construct matrix D
    d = Dic(:,sup(i));
    D = [D d];
    
%     for j=1:i
%        D(:,j) = D(:,j)./sum(abs(D(:,j))); 
%     end


    % 4) amplitude adjustment (least squares)
    a = pinv(D)*signal;
    %a = ((D.'*D)^(-1))*D.'*signal;
    %a = D.'*signal;

    % 5) new residual calculation
    ajuste(:,i) = D*a;
    %r(:,i+1) = signal - ajuste(:,i);
    r = signal - ajuste(:,i);

end

%     figure
%     subplot(3,3,1), hold on, grid on, stem(r(:,1)), title('r: 0')
%     subplot(3,3,2), hold on, grid on, stem(abs(r(:,1))), title('abs(r): 0')
%     subplot(3,3,3), hold on, grid on, stem(angle(r(:,1))), title('angle(r): 0')
%     subplot(3,3,4), hold on, grid on, stem(r(:,2)), title('r: 1')
%     subplot(3,3,5), hold on, grid on, stem(abs(r(:,2))), title('abs(r): 1')
%     subplot(3,3,6), hold on, grid on, stem(angle(r(:,2))), title('angle(r): 1')
%     subplot(3,3,7), hold on, grid on, stem(r(:,3)), title('r: 2')
%     subplot(3,3,8), hold on, grid on, stem(abs(r(:,3))), title('abs(r): 2')
%     subplot(3,3,9), hold on, grid on, stem(angle(r(:,3))), title('angle(r): 2')
% 
%     figure                                                 
%     subplot(3,2,1), hold on, grid on, stem(Dic(:,sup(1)-1)), title(['Atomo: ' num2str(sup(1)-1)]);
%     subplot(3,2,2), hold on, grid on, stem(Dic(:,sup(2)-1)), title(['Atomo: ' num2str(sup(2)-1)]);
%     subplot(3,2,3), hold on, grid on, stem(Dic(:,sup(1))), title(['Atomo: ' num2str(sup(1))]);
%     subplot(3,2,4), hold on, grid on, stem(Dic(:,sup(2))), title(['Atomo: ' num2str(sup(2))]);
%     subplot(3,2,5), hold on, grid on, stem(Dic(:,sup(1)+1)), title(['Atomo: ' num2str(sup(1)+1)]);
%     subplot(3,2,6), hold on, grid on, stem(Dic(:,sup(2)+1)), title(['Atomo: ' num2str(sup(2)+1)]);
% 
% %     figure                                                 
% %     subplot(3,2,1), hold on, grid on, stem(abs(Dic(:,sup(1)-1))), title(['Atomo: ' num2str(sup(1)-1)]);
% %     subplot(3,2,2), hold on, grid on, stem(abs(Dic(:,sup(2)-1))), title(['Atomo: ' num2str(sup(2)-1)]);
% %     subplot(3,2,3), hold on, grid on, stem(abs(Dic(:,sup(1)))), title(['Atomo: ' num2str(sup(1))]);
% %     subplot(3,2,4), hold on, grid on, stem(abs(Dic(:,sup(2)))), title(['Atomo: ' num2str(sup(2))]);
% %     subplot(3,2,5), hold on, grid on, stem(abs(Dic(:,sup(1)+1))), title(['Atomo: ' num2str(sup(1)+1)]);
% %     subplot(3,2,6), hold on, grid on, stem(abs(Dic(:,sup(2)+1))), title(['Atomo: ' num2str(sup(2)+1)]);
% %  
% %     figure                                                 
% %     subplot(3,2,1), hold on, grid on, stem(angle(Dic(:,sup(1)-1))), title(['Atomo: ' num2str(sup(1)-1)]);
% %     subplot(3,2,2), hold on, grid on, stem(angle(Dic(:,sup(2)-1))), title(['Atomo: ' num2str(sup(2)-1)]);
% %     subplot(3,2,3), hold on, grid on, stem(angle(Dic(:,sup(1)))), title(['Atomo: ' num2str(sup(1))]);
% %     subplot(3,2,4), hold on, grid on, stem(angle(Dic(:,sup(2)))), title(['Atomo: ' num2str(sup(2))]);
% %     subplot(3,2,5), hold on, grid on, stem(angle(Dic(:,sup(1)+1))), title(['Atomo: ' num2str(sup(1)+1)]);
% %     subplot(3,2,6), hold on, grid on, stem(angle(Dic(:,sup(2)+1))), title(['Atomo: ' num2str(sup(2)+1)]);
%     
%     figure  
%     subplot(2,3,1), hold on, grid on, stem(ajuste(:,1)), title('D*a')
%     subplot(2,3,2), hold on, grid on, stem(abs(ajuste(:,1))), title('abs(D*a)')
%     subplot(2,3,3), hold on, grid on, stem(angle(ajuste(:,1))), title('anle(D*a)')
%     subplot(2,3,4), hold on, grid on, stem(ajuste(:,2))
%     subplot(2,3,5), hold on, grid on, stem(abs(ajuste(:,2)))
%     subplot(2,3,6), hold on, grid on, stem(angle(ajuste(:,2)))
    
for i=1:length(a)
    amps(sup(i)) = a(i);
    
end


% sup = sup;

end