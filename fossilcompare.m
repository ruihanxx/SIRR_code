cla
clc;
clear
m = 5000;
n = 200;
s = 600;
conds = [-1:-1:-15];
res_sizes = -conds-16;
F1f = zeros(1,15);
F2f = zeros(1,15);
F1b = zeros(1,15);
F2b = zeros(1,15);
for j=1:length(conds)
    k = conds(j);
    e_length = 10^(res_sizes(j));
    %generate matrix
    R = normrnd(0,1,m,n);
    [U,~] = qr(R,0);
    R = normrnd(0,1,n,n);
    [V,~] = qr(R,0);
    
    Sigma = diag(logspace(0,k,n));
    
    A = U*Sigma*transpose(V);
    x0 = normrnd(0,1,n,1);
    e = normrnd(0,1,m,1);
    e = e-U*transpose(U)*e;
    e = e/norm(e)*e_length;
    b = A*x0 + e;
    xstar = x0;
    %solve
    t1 = 9;
    J = 3;
    t2 = max(round(log(k/log(0.5))/log(J))+1,1);
    J1 = 8;
    
    K = 3;
    
    S = normrnd(0,1,s,m)/sqrt(s);
    
    [~,x_list1,time_list1] = IDS_solver(A,b,s,K,[J,J1],t2,1e-17);
    [~,x_list2,~,~] = fossils(A,b,s,[],[],true);
    x1 = x_list1(:,end);
    x2 = x_list2(:,end);
    
    F1f(j) = [norm(x1-xstar)/norm(xstar)];
    F2f(j) = [norm(x2-xstar)/norm(xstar)];
    
    F1b(j) = [backward_error_ls(A,b,x1)];
    F2b(j) = [backward_error_ls(A,b,x2)];
end




if true
figure(1);
F1 = F1f;
F2 = F2f;
loglog(10.^(-conds),F1,'--*','LineWidth', 3)
hold on
loglog(10.^(-conds),F2,'--o','LineWidth', 3)
hold on


legend({'SIRR','FOSSILS'},'Location','best','FontSize',20)
xlim([10^1 10^15])


ylabel('forward error (log_{10})','FontSize',20)

xlabel('$difficulty=\kappa=\frac{\|b-Ax^*\|}{u}$ ','interpreter','latex','FontSize',20)

base = '/Users/carpeduem/Desktop/learning/deep learning/project/randomness in algebra/my/simulation/paper/picture/fossil compare/';
fn = strcat(base,'fossil_fe.png');
saveas(figure(1),fn)
end

if true
figure(2);
F1 = F1b;
F2 = F2b;
loglog(10.^(-conds),F1,'--*','LineWidth', 3)
hold on
loglog(10.^(-conds),F2,'--o','LineWidth', 3)
hold on


legend({'SIRR','FOSSILS'},'Location','best','FontSize',20)
xlim([10^1 10^15])


ylabel('backward error (log_{10})','FontSize',20)

xlabel('$difficulty=\kappa=\frac{\|b-Ax^*\|}{u}$ ','interpreter','latex','FontSize',20)

base = '/Users/carpeduem/Desktop/learning/deep learning/project/randomness in algebra/my/simulation/paper/picture/fossil compare/';
fn = strcat(base,'fossil_be.png');
saveas(figure(2),fn)
end
