% CORRELAÇÃO DE AMPLITUDE
% DICIONÁRIO DE AMPLITUDES


function [amps, a, sup, res] = CORR_amp(k, signal, Dic);

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
%      temp = Dic(3:255,:)'*r(3:255,i);
%      [~, index] = max(abs(temp));
     
    center_r = sum(real(r(3:255)));
    
    %figure, hold on, grid on, stem(real(r(1:255)));
    %figure, hold on, grid on, stem(sum(real(Dic(3:255,)));
    
    if center_r > 0
        [rho] = corr(abs(r),abs(Dic(:,1:512)),'type', 'Pearson'  );
        [~, index] = max(rho(1:end).^10);
    elseif center_r <0
        [rho] = corr(abs(r),abs(Dic(:,512:end)),'type', 'Pearson'  );
        [~, index] = max(rho(1:end).^10);
        index = index + 512 -1;
    end
    
%    [rho,pval] = corr(abs(r(2:256)),abs(Dic(2:256,:)),'type', 'Pearson'  );
    %figure, stem(rho), title('abs(rho)');
    %figure, stem(pval), title('abs(pval)'); 
    
%    [~, index] = max(rho(1:end).^2);
    
    
     sup(i) = index(1);

    % 3) construct matrix D
    d = Dic(:,sup(i));
    D = [D d];
    
%    D = [d conj(d)];
%     if i==1
%         D = [D d]
%     else
%         D = [D conj(d));
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
for i=1:1;
    amps(sup(i)) = a(i);
    
end

res = r;
% sup = sup;

end