%%%% OMP
% implementação propria, ajustado para dicionário freq domain


function [amps, sup, res] = OMP_freq4(k, signal, Dic)

[m,n] = size(Dic); % Get dictionary dimension

% 1) initializations
error_goal = 1e-3;
amps = zeros(1,n);
%a = zeros(1,m);
% sup = zeros(1,m);
D = [];

r = signal;   % residue
%norm = sum(abs(r));
%r = r./norm;
% r = r.^2;
% Dic = Dic.^2;

i=0;

%while (error > error_goal) || (i<20),
while (i<k)
    
    i = i+1;
    % 2) obtain support vector (search indexes) 
    isodd = mod(i,2);
    
    if isodd==1
        %temp = Dic(1:m/2,:)'*r(1:m/2,i);
        temp = Dic'*r(:,i);
        [~, index] = max(abs(temp(1:n/2)));
        sup(i) = index(1);
        d = Dic(:,sup(i));
        D = [D d];
    else
        %temp = Dic(m/2+1:end,:)'*r(m/2+1:end,i);
        %[~, index] = max(abs(temp(n/2+1:end)));
        %sup(i) = index(1)+n/2;
        sup(i) = n-sup(i-1)+1;
        %d = Dic(:,sup(i));
        %D = [D d];
        D = [D conj(flip(d))];
    end
%     figure
%     subplot(2,1,1), grid on, hold on, stem(real(temp)), title('real(temp)')
%     subplot(2,1,2), grid on, hold on, stem(imag(temp)), title('imag(temp)')
%     figure
%     subplot(2,1,1), grid on, hold on, stem(abs(temp)), title('abs(temp)')
%     subplot(2,1,2), grid on, hold on, stem(angle(temp)), title('angle(temp)')

    
    

%     % 3) construct matrix D
%     d = Dic(:,sup(i));
%     D = [D d];
    
    %d1 = Dic(:,sup(i));
    %d2 = Dic(:,n-sup(i)+1);
    %D = [D d1 d2];
% % 1) Dois vetores:
%    D = [d conj(flip(d))];
% % 2) concatenar em um único vetor (Nx1) com seu próprio conjugado:   
%    D = [d(1:m/2+1); conj(flip(d(2:m/2)))];
% % 3) adicionando a imagem extraída do dicionário:
%     d1 = Dic(:,30);
%     d2 = Dic(:,2788);
%     D = [d1 d2];
%     d1 = Dic(1:m/2+1,29);
%     d2 = Dic(m/2+2:end,2789);
%     d1 = Dic(4,29);
%     d2 = Dic(254,2789);
%     D = [d1 d2];
%     a = D'*[signal(4) signal(254)];

    % 4) amplitude adjustment (least squares)
    a = pinv(D)*signal;
    %a = ((D'*D)^(-1))*D'*signal;
    %a = D'*signal;
%     a = pinv(D(1:m/2,:))*signal(1:m/2);
    
    %a = sum(D'*signal)/sum(D.^2);
%     a(1,1) = signal(4)/D(4,1);
%     a(2,1) = signal(254)/D(254,2);

    % 5) new residual calculation
     ajuste(:,i) = D*a;
     r(:,i+1) = signal - ajuste(:,i);
%     norm = sum(abs(r(:,i+1)));
%     r(:,i+1) = r(:,i+1)./norm;

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
%     
for i=1:k;
    amps(sup(i)) = a(i);
    
end

res = r;
% sup = sup;

end