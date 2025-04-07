
M = csvread('../data/susy.csv',1,0);
X = M(1:1e6,2:end);
b = M(1:1e6,1);
bandwidth = 4;
S = randsample(1e6,1e3);
A = exp(-pdist2(X,X(S,:),"euclidean").^2 / (2*bandwidth^2));

sizes = round(logspace(1,3,11));
Time = zeros(3,11);
E = zeros(3,11);
E2 = zeros(3,11);
Cond = [];
for i = 1:11
    sizes(i);
    AA = A(:,1:sizes(i));
    n = zeros(1,500);
    cond = sqrt(max(eig(AA'*AA)))
    Cond(end+1) = cond;
    % QR
    tic; xqr = AA\b; Time(1,i) = toc;
    % E2(1,i) = backward_error_ls(AA, b, xqr);
    %fossils
    tic; xf = fossils(AA,b,6*sizes(i)); Time(2,i) = toc;
    % E2(2,i) = backward_error_ls(AA, b, xf);
    %SIRR
    J1 = 20;
    J = 3;
    t = round(log(log(1/cond)/log(0.25))/log(J));
    tic; xs = IDS_solver_svd(AA,b,6*sizes(i),3,[J,J1],t,5e-16); Time(3,i) = toc;
    % E2(3,i) = backward_error_ls(AA, b, xs);
    
end
