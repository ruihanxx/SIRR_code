cla
clc;
clear
m = 2000;
n = 50;
s = 200;
conds = [-8];
res_sizes = [-4,-3,-2,-1,-0.5,0,0.5,1,1.5,2];
ratio = [];
k = conds(1);


F1b = [];
F2b = [];
F3b = [];

F1f = [];
F2f = [];
F3f = [];
for j=1:length(res_sizes)



%generate matrix
R = normrnd(0,1,m,n);
[U,~] = qr(R,0);
R = normrnd(0,1,n,n);
[V,~] = qr(R,0);

Sigma = diag(logspace(0,k,n));

A = U*Sigma*transpose(V);
x0 = normrnd(0,1,n,1);
e_length = 10^(res_sizes(j))*norm(A*x0);
e = normrnd(0,1,m,1);
e = e-U*transpose(U)*e;
e = e/norm(e)*e_length;
b = A*x0 + e;

xstar = x0;
lambda = 0;
ratio(end+1) = e_length/norm(A*xstar);

method1_result = cell(length(conds), length(res_sizes));
method2_result = cell(length(conds), length(res_sizes));
qr_result = cell(length(conds), length(res_sizes));

%xstar = (transpose(A)*A+lambda*eye(n))\(transpose(A)*b);
%solve
t1 = 9;
t2 = max(round(log2(k/log(0.5)))+1,1);
J = 3;
J1 = 5;
K = 4;
num = 1;

E1 = [];
E2 = [];
E3 = [];

S = normrnd(0,1,s,m)/sqrt(s);

[x1,x_list1,time_list1] = IDS_solver(A,b,s,K,[J,J1],t2);
x0=x_list1(:,2)*0;
x = IDS_solver_test(A,b,S,5,3,t1,x0);



xqr = A\b;
F1f(end+1) = [norm(x1-xstar)/norm(xstar)];
F2f(end+1) = [norm(x-xstar)/norm(xstar)];
F3f(end+1) = [norm(xqr-xstar)/norm(xstar)];

F1b(end+1) = [backward_error_ls(A,b,x1)];
F2b(end+1) = [backward_error_ls(A,b,x)];
F3b(end+1) = [backward_error_ls(A,b,xqr)];
end




if true
figure(1);
F1 = F1f;
F2 = F2f;
F3 = F3f;
plot(log10(ratio),log10(F1),'LineWidth', 3)
hold on
plot(log10(ratio),log10(F2),'LineWidth', 3)
hold on
plot(log10(ratio),log10(F3),'.--','LineWidth', 3)
hold on
% yline(-16,'-','machine precision')
yl = [-17,10];
patch([0, max(log10(ratio)),max(log10(ratio)), 0], [yl(1), yl(1),yl(2),yl(2)], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
% plot(log(E3(2,:)),'LineWidth', 3)
legend({'SIRR','SRR','mldivide(matlab)'},'Location','best','FontSize',20)
xlim([min(log10(ratio)) max(log10(ratio))])
ylim(yl)

ylabel('forward error (log_{10})','FontSize',20)
title(strcat('condition number = 1e',num2str(-k)),'FontSize',20)
xlabel('$log_{10}(\frac{\|b-Ax^*\|}{\|Ax^*\|})$ ','interpreter','latex','FontSize',20)

base = '/Users/carpeduem/Desktop/learning/deep learning/project/randomness in algebra/my/simulation/paper/picture/error_size/';
fn = strcat(base,'fe',num2str(-k),'.png');
saveas(figure(1),fn)
end

if true
figure(2);
F1 = F1b;
F2 = F2b;
F3 = F3b;
plot(log10(ratio),log10(F1),'LineWidth', 3)
hold on
plot(log10(ratio),log10(F2),'LineWidth', 3)
hold on
plot(log10(ratio),log10(F3),'.--','LineWidth', 3)
hold on
% yline(-16,'-','machine precision')
yl = [-17,-10];
patch([0, max(log10(ratio)),max(log10(ratio)), 0], [yl(1), yl(1),yl(2),yl(2)], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
% plot(log(E3(2,:)),'LineWidth', 3)
legend({'SIRR','SRR','mldivide(matlab)'},'Location','best','FontSize',20)
xlim([min(log10(ratio)) max(log10(ratio))])
ylim(yl)

ylabel('backward error (log_{10})','FontSize',20)
title(strcat('condition number = 1e',num2str(-k)),'FontSize',20)
xlabel('$log_{10}(\frac{\|b-Ax^*\|}{\|Ax^*\|})$ ','interpreter','latex','FontSize',20)

base = '/Users/carpeduem/Desktop/learning/deep learning/project/randomness in algebra/my/simulation/paper/picture/error_size/';
fn = strcat(base,'be',num2str(-k),'.png');
saveas(figure(2),fn)
end