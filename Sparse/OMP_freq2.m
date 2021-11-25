%%%% OMP
% implementação propria, ajustado para dicionário freq domain

% MODIFICAÇÃO NA SELEÇÃO DO ÁTOMO: TESTE DE MSE NO RESÍDUO

function [amps, sup, r] = OMP_freq2(k, signal, Dic, GranH, P, freq);

[m,n] = size(Dic); % Get dictionary dimension

 % AVALIAR CRITÉRIO PARA DETERMINAÇÃO DA MARGEM

% 1) initializations
error_goal = 1e-3;
amps = zeros(1,n);
%a = zeros(1,m);
sup = zeros(1,k);       % vetor de suporte
% % margin = floor(P/2);    % determina margem de índices admissível
% r1  = ones(m,margin);   % vetor de residuos
% a_mse_r1 = ones(1,margin);% vetor do MSE para o resíduo

% for j=1:n
%     gama(j) = sum(abs(Dic(:,j)));      
% end
% Inicializações
r = signal;             % residue
% fator = sum(signal);
% norm_r = signal./fator;
D = [];

%     for l=1:n 
%         %soma(l) = sum(abs(Dic(:,l)));
%         %Dic_n(:,l) = Dic(:,l)./soma(l);
%         Dic_n(:,l) = Dic(:,l)./abs(Dic(:,l));
%     end

i=0;
% Delta_f = freq - f0;
%while (error > error_goal) || (i<20),
while (i<k)
    
    i = i+1;
    % 2) obtain support vector (search indexes)
    
    a = Dic'*r;
    [~, index] = max(abs(a));
    ac = abs(a).*cos(angle(a));
    as = abs(a).*sin(angle(a));
    [proj, index_c] = max(abs(ac));
    [erro, index_s] = min(abs(as));
    sup(i) = index(1);
    
    figure, stem(angle(a)), title('proj');
     figure, stem(as), title('erro');

%     if i==1
%         sup(i) = 29; %P+fix(Delta_f/GranH)+1;
%         %sup(i) = P+1;%fix(Delta_f*2)+1;
%     else
        %sup(2) = 619;
%         %temp = Dic.'*r(:,i);
%         a = Dic.'*r;
%         [~, index] = max(abs(a));
%         aa = abs(a).*cos(angle(a));
%         [proj, index2] = max(abs(aa));
%         sup(i) = index(1);
%    end
    
%     figure
%     subplot(3,1,1), grid on, hold on, stem((a)), title('a = Dic.''*r')
%     subplot(3,1,1), grid on, hold on, stem(index,(a(index)))
%     subplot(3,1,2), grid on, hold on, stem(abs(a)), title('abs(a)')
%     subplot(3,1,2), grid on, hold on, stem(index,abs(a(index)))
%     subplot(3,1,3), grid on, hold on, stem(angle(a)), title('angle')
%     subplot(3,1,3), grid on, hold on, stem(index,angle(a(index)))
%    
%     figure
%     subplot(3,1,1), grid on, hold on, stem(ac), title('abs(a)*cos(angle(a))')
%     subplot(3,1,1), grid on, hold on, stem(index_c,proj)
%     subplot(3,1,2), grid on, hold on, stem(abs(ac)), title('abs(abs(a)*cos(angle(a)))')
%     subplot(3,1,2), grid on, hold on, stem(index_c,abs(proj))
%     subplot(3,1,3), grid on, hold on, stem(angle(ac)), title('angle(abs(a)*cos(angle(a)))')
%     subplot(3,1,3), grid on, hold on, stem(index_c,angle(proj))
%     
%     figure
%     subplot(3,1,1), grid on, hold on, stem(as), title('abs(a)*sin(angle(a))')
%     subplot(3,1,1), grid on, hold on, stem(index_s,erro)
%     subplot(3,1,2), grid on, hold on, stem(abs(as)), title('abs(abs(a)*sin(angle(a)))')
%     subplot(3,1,2), grid on, hold on, stem(index_s,abs(erro))
%     subplot(3,1,3), grid on, hold on, stem(angle(as)), title('angle(abs(a)*sin(angle(a)))')
%     subplot(3,1,3), grid on, hold on, stem(index_s,angle(erro))
%     

    
    
    
%     % Itera -margin/2 ate margin/2 átomos do escolhido anteriormente (index);
%     lim = floor(margin/2)+1;             % define limite para inicio da iteração
%     %figure
% 
%     for j=1:margin;
% 
%         index1(j) = sup(i)-lim+j;
%         r1(:,j) = r - Dic(:,index1(j));   % calcula residuo para cada átomo
%         abs_r1(:,j) = abs(r) - abs(Dic(:,index1(j)));   % calcula residuo para cada átomo
%         ang_r1(:,j) = angle(r) - angle(Dic(:,index1(j)));   % calcula residuo para cada átomo
%         
%         for n=1:m % corrige angulos
%             if ang_r1(n,j) > pi;
%                 ang_r1(n,j) = ang_r1(n,j) - pi;
%             elseif ang_r1(n,j) < -pi;
%                 ang_r1(n,j) = ang_r1(n,j) + pi;
%             end
%         end
%         %ang_r1(:,j) =ang_r1(:,j)./pi;
%         tve(:,j) =  (abs(r1(:,j)))./abs(Dic(:,index1(j)));
%         
%         mse_tve(:,j) = mean(abs(tve(:,j)));
%         
%         mse_r1(:,j) =  mean(abs(r1(:,j)).^2); 
%         mse_abs_r1(:,j) =  mean(abs(abs_r1(:,j)).^2);
%         mse_ang_r1(:,j) =  mean(abs(ang_r1(:,j)).^2); 
%         
%         m_ang_r1(:,j) =  mean(angle(r1(:,j)).^2);% calcula MSE de cada resíduo
%         c_mse_r1(:,j) =  mean((abs(r1(:,j)).*angle(r1(:,j))).^2);    % calcula MSE de cada resíduo
% 
%         %subplot(5,1,j), hold on, grid on, stem(Dic(:,index1(j))), title(['indice: ', num2str(index1(j))])
%     end
%     

% %     
% %     figure
% %     subplot(3,1,1), grid on, hold on, stem(mse_r1), title(['mean(abs(r).^2) K = ', num2str(i)])
% %     subplot(3,1,2), grid on, hold on, stem(m_ang_r1), title(['mean(angle)^2 - K = ', num2str(i)])
% %     subplot(3,1,3), grid on, hold on, stem(c_mse_r1), title(['mean(abs*angle)^2 - K = ', num2str(i)])
% %     
% %     figure
% %     subplot(3,1,1), grid on, hold on, stem(mse_tve),  title(['mean(tve) K = ', num2str(i)])
% %     subplot(3,1,2), grid on, hold on, stem(mse_abs_r1),  title(['mean(abs-abs) K = ', num2str(i)])
% %     subplot(3,1,3), grid on, hold on, stem(mse_ang_r1),  title(['mean(ang-ang) K = ', num2str(i)])    
% 
%     [~, index] = min(mse_r1);
%     sup(i) = index1(1)+ index(1)-1;
%     figure
%     subplot(3,1,1), hold on, grid on
%     stem(r(:,i))
%     subplot(3,1,2), hold on, grid on
%     stem(Dic(:,10))
%     subplot(3,1,3), hold on, grid on
%     stem(Dic(:,12))
%     
%     figure
%     plot(abs(temp))
    
    % 3) construct matrix D
    d = Dic(:,sup(i));
    D = [D d];
    
%     for j=1:i
%        D(:,j) = D(:,j)./sum(abs(D(:,j))); 
%     end


    % 4) amplitude adjustment (least squares)
    amp = pinv(D)*signal;
    %a = ((D.'*D)^(-1))*D.'*signal;
    %a = D.'*signal;

    % 5) new residual calculation
    ajuste(:,i) = D*amp;
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
    
for i=1:length(amp)
    amps(sup(i)) = amp(i);
    
end


% sup = sup;

end