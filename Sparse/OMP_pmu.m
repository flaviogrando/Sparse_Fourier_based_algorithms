%%%% OMP



function [amps , sup, res] = OMP_pmu(k, signal, Dic);

% 1) initializations
r = signal;   % residue
D= [];

% 2) obtain support vector (search indexes)
for i=1:k
    [v , index] = max(Dic.'*r);
    sup(i) = index(1);
end

% 3) construct matrix D
for i=1:length(sup(i))
    
    d = Dic(:,sup(i));
    D = [D d];
end

% 4) amplitude adjustment (least squares)
a = pinv(D)*signal;


% 5) new residual calculation
res = signal - D*a;
amps = a;

end