% CORRELAÇÃO DE AMPLITUDE
% DICIONÁRIO DE AMPLITUDES

%%%% OMP
% implementação propria, ajustado para dicionário freq domain


function [amps, a, sup, res] = CORR_amp2(k, signal, Dic);

%load('Dic_A.mat')

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

%while (error > error_goal) || (i<20),
while (i<k)
    
    i = i+1;
    % 2) obtain support vector (search indexes)
    %temp = Dic'*r(:,i);
    %temp = temp/abs(r(:,i));
    
    center_r = sum(real(r(3:255)));
    
    
    
    [rho,pval] = corr(real(r(2:256,i)),real(Dic(2:256,:)),'type', 'Pearson'  );

    for j=1:n
        r = r/sum(abs(r));
        erro(:,j) = (Dic(1:256,j) - r(1:256))./sum(abs(Dic(1:256,j)));
        sum_erro(j) = sum(abs(erro(:,j)));
    end
%     figure, mesh(abs(erro)), title('erro');
%     figure, stem(sum_erro)
    
    [~, index] = min(sum_erro);
    sup(i) = index(1);
%     figure, stem(rho), title('abs(rho)');
%     figure, stem(pval), title('abs(pval)'); 
    
%     if i==1
% %         temp = Dic(1:m/2,:)'*r(1:m/2,i);
% %         [~, index] = max(abs(temp(1:n/2)));
%         
%         [~, index] = max(rho(1:end).^2);
%          %index = index+512-1;
%         sup(i) = index(1);
%     else
% %         temp = Dic(m/2+1:end,:)'*r(m/2+1:end,i);
% %         [~, index] = max(abs(temp(n/2+1:end)));
%         %[rho,pval] = corr(abs(r),abs(Dic(:,:)),'type', 'Pearson'  );
%         [~, index] = max(rho(257:end));
%         sup(i) = index(1);%+n/2;  
%     end
    

    
%     figure
%     subplot(2,1,1), grid on, hold on, stem(real(temp(1:50))), title('real(temp)')
%     subplot(2,1,2), grid on, hold on, stem(imag(temp(1:50))), title('imag(temp)')
%     figure
%     subplot(2,1,1), grid on, hold on, stem(abs(temp(1:50))), title('abs(temp)')
%     subplot(2,1,2), grid on, hold on, stem(angle(temp(1:50))), title('angle(temp)')

%     [rho,pval] = corr(abs(r),abs(Dic(:,:)),'type', 'Spearman'  );
%     figure, stem(rho), title('abs(rho)');
%     figure, stem(pval), title('abs(pval)');
%     
%     [rho,pval] = corr(angle(r),angle(Dic(:,:)));
%     figure, stem(rho), title('angle(rho)');
%     figure, stem(pval), title('angle(pval)');
    
%     [rho,pval] = corr(abs(r),abs(Dic(:,index(1)-2:index(1)+2)),'type', 'Pearson');
%     %[rho,pval] = corr(angle(r),angle(Dic(:,index(1)-2:index(1)+2)));
% 
%     [~, index2] = min(rho);
%     [~, index3] = max(pval);

%      index2 = index2+index(1)-3;
%      index3 = index3+index(1)-3;
    
%     sup(i) = index(1);
% %     
%     figure, stem(rho), title('abs(rho)');
%     figure, stem(pval), title('abs(pval)');
%     if i==1
%         sup(i) = 29; %round(freq/gran) + 2;
%     end
    
    % 3) construct matrix D
    d = Dic(:,sup(i));
    D = [D d];
    
%    D = [d conj(d)];
%     if i==1
%         D = [D d]
%     else
%         D = [D conj(d));
%     end
    
    
%     for j=1:i
%        D(:,j) = D(:,j)./sum(abs(D(:,j))); 
%     end


    % 4) amplitude adjustment (least squares)
    a = pinv(D)*signal;
    %a = ((D.'*D)^(-1))*D.'*signal;
    %a = D.'*signal;

    % 5) new residual calculation
    ajuste(:,i) = D*a;
    r(:,i+1) = r(:,i) - ajuste(:,i);
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