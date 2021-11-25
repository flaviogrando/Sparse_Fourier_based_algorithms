%%%% OMP
% implementação propria, ajustado para dicionário freq domain


function [amps, sup, res] = OMP_freq1(k, signal, Dic)

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
% r = r.^2;
% Dic = Dic.^2;

i=0;
erro = sum(abs(r));

%while erro > 0.1;
while (i<k)
    
    i = i+1;
    % 2) obtain support vector (search indexes)
%     temp = Dic'*r(:,i);
%     [~, index] = max(abs(temp));
%     sup(i) = index(1);
    
    isodd = mod(i,2);
  
    if isodd==1
        temp = Dic(1:m/2,1:n/2)'*r(1:m/2,i);
        [~, index] = max(abs(temp(1:n/2)));
        sup(i) = index(1);
%         temp = Dic'*r(:,i);
%         [~, index] = max(abs(temp(1:n/2)));
%         sup(i) = index(1);
    else
        
%         temp = Dic(m/2+1:end,:)'*r(m/2+1:end,i);
%         [~, index] = max(abs(temp(n/2+1:end)));
%         sup(i) = index(1)+n/2;
        temp = Dic(1:end,n/2:end)'*r(1:end,i);
        [~, index] = max(abs(temp));
        sup(i) = index(1)+n/2;
    end
    
%     lim = n/2;
%     figure
%     subplot(2,1,1), grid on, hold on, stem(real(temp(1:lim))), title('real(temp)')
%     subplot(2,1,2), grid on, hold on, stem(imag(temp(1:lim))), title('imag(temp)')
%     figure
%     subplot(2,1,1), grid on, hold on, stem(abs(temp(1:lim))), title('abs(temp)')
%     subplot(2,1,2), grid on, hold on, stem(angle(temp(1:lim))), title('angle(temp)')
% % % 
%     if isodd==1
%         figure,grid on, hold on,
%         stem(abs(temp(1:end)))
%         xlabel('índice'), ylabel('amplitude (p.u.)')
%     else
%         stem(abs(flip(temp(1:end))))
%     end

    
    % 3) construct matrix D
    d = Dic(:,sup(i));
    D = [D d];
    

    % 4) amplitude adjustment (least squares)
    a = pinv(D)*signal;
    %a = ((D'*D)^(-1))*D'*signal;
    %a = D'*signal;
%     if isodd==1
%         a = pinv(D(1:m/2,:))*signal(1:m/2);
%     else
%         a = pinv(D)*signal;
%     end

    % 5) new residual calculation
    ajuste(:,i) = D*a;
    r(:,i+1) = signal - ajuste(:,i);
    %r(:,i+1) = r(:,i) - ajuste(:,i);
    
    %erro = sum(abs(r(:,i+1)));
%     norm = sum(abs(r(:,i+1)));
%     r(:,i+1) = r(:,i+1)./norm;
% if i>=10-1
%     erro = 0;
% end

%%% Extra) ajuste fino


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
     
for i=1:k;
    
    amps(sup(i)) = a(i);
        
end

res = r;
% sup = sup;

end