cla
clc;
clear
m = 2000;
n = 100;
S = [300,350,400,450,500,550,600,650,700,750,800]/2;
conds = [-12];
res_sizes = [-3];
k = conds(1);

e_length = 10^(res_sizes(1));

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
lambda = 0;
F1 = [];
F2 = [];
for j=1:length(S)

s=S(j);

method1_result = cell(length(conds), length(res_sizes));
method2_result = cell(length(conds), length(res_sizes));
qr_result = cell(length(conds), length(res_sizes));

%xstar = (transpose(A)*A+lambda*eye(n))\(transpose(A)*b);
%solve
t1 = 10;
t2 = round(log2(k/log(0.5)))+1;
J = 2;
J1 = 5;
K = 3;
num = 100;

E1 = [];
E2 = [];
E3 = [];
for i=1:num
x_list1 = inf;





[x1,x_list1,time_list1] = IDS_solver(A,b,s,K,[J,5],t2);
[x2,x_list2,~,~] = fossils(A,b,s,[],[],true);

eA1t = norm(A'*(A*x1-b))/norm(A'*b);
eA2t = norm(A'*(A*x2-b))/norm(A'*b);



E1(end+1) = [eA1t];
E2(end+1) = [eA2t];


end
F1(end+1)=sum(E1>1e-5)/num;
F2(end+1)=sum(E2>1e-5)/num;
end
% ans1=sum(E1,2)/num
% ans2=sum(E2,2)/num

plot(S(1:end)/n,F1(1:end),'LineWidth', 3)
hold on
plot(S(1:end)/n,F2(1:end),'LineWidth', 3)
hold on
% plot(log(E3(2,:)),'LineWidth', 3)
legend({'SIRR','FOSSILS'},'FontSize',20)

ylabel('fail rate','FontSize',20)
xlabel('sketch dimension:s/n','FontSize',20)