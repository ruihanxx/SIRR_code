cla
clc;
clear
M = [2000,4000,6000,8000,10000];
Num = zeros(3,0);
Time = zeros(3,0);
for j=1:length(M)
m = M(j);
n = M(j)*0.05;
s = M(j)*0.2;
lambda = 0;

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
%xstar = (transpose(A)*A+lambda*eye(n))\(transpose(A)*b);
%solve
t1 = 10;
t2 = round(log2(max(k,-16-k)/log(0.5)));
J = 2;
J1 = 10;
K = 4;
num = 30;
num1 = num;
num2 = num;
num3=num;
time_avg1 = 0;
time_avg2 = 0;
time_avg3 = 0;


for i = 1:num
tic
[x,x_list1] = IDS_solver_svd(A,b,s,K,[J,J1],t2,5e-15);
time1=toc;
% if backward_error_ls(A,b,x_list1(:,end))<1e-15
    time_avg1 = time_avg1+time1;
% else
%     num1=num1-1;
% end

tic
[~,x_list2,~,~] = fossils(A,b,s,[],[],true);
time2=toc;
% if backward_error_ls(A,b,x_list2(:,end))<1e-15
    time_avg2 = time_avg2+time2;
% else
%     num2=num2-1;
% end

tic
A\b;
time3=toc;
% if backward_error_ls(A,b,x_list2(:,end))<1e-15
    time_avg3 = time_avg3+time3;
% else
%     num2=num2-1;
% end
end
time_avg1 = time_avg1/num1;
time_avg2 = time_avg2/num2;
time_avg3 = time_avg3/num3;
Time(:,end+1) = [time_avg1;time_avg2;time_avg3]
Num(:,end+1) = [num1;num2;num3];
end

plot(M,Time(1,:),'LineWidth',2)
hold on
plot(M,Time(2,:),'LineWidth',2)
hold on 

legend({'SIRR','FOSSILS','mldivide'},'Location','best','FontSize',20)
ylabel('runtime','FontSize',14,'FontSize',20)
xlabel('demension of A','FontSize',20)