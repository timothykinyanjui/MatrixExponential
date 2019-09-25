% Script to integrate the full system - plots simple infected curve rather
% than complex figure in main paper

clear

params.beta = 0.5;
params.N = 100;
params.g = 1;
params.mu = 0;
params.alpha = 0.06;
tmax = 10;
i0 = 1e-3;

[H,ss,ii] = SIRS(params.N,params.beta/(params.N-1),params.g,params.mu,0,[0 params.N],[0 params.N]);
K = SIRS(params.N,0,0,0,params.alpha,[0 params.N],[0 params.N]);

lp = length(ii);
params.ii = ii;
pinit = zeros(length(ii),1);
pinit(ii==1 & ss==(params.N-1))=i0;
pinit(ii==0 & ss==(params.N))=1-i0;

Qdiff = @(t,p) (H + (((ii')*p))*K)*p;

h = 0.1;
tran = 0:h:tmax;
lt = length(tran);

tic
[T, P] = ode45(Qdiff, tran, pinit);
t45 = toc;

%%

I = sparse(diag(ones(lp,1)));

p1 = zeros(lp,lt);
tic
p1(:,1) = pinit;
for i=2:lt
    M = (H + (((ii')*p1(:,i-1)))*K);
    A1 = I - h*M;
    p1(:,i) = A1\(p1(:,i-1));
end
t1 = toc;
P1 = p1';

p2 = zeros(lp,lt);
tic
p2(:,1) = pinit;
for i=2:lt
    M = (H + (((ii')*p2(:,i-1)))*K);
    A1 = I - h*M + (1/2)*(h^2)*(M^2);
    ptemp = A1\(p2(:,i-1));
    Mdot = (((ii')*M*p2(:,i-1)))*K;
    A2 = I - (1/2)*(h^2)*Mdot;
    p2(:,i) = A2\(ptemp);
end
t2 = toc;
P2 = p2';
% 
% A3 = sparse(diag(ones(lp,1)) - h*M + ...
%     (1/2)*(h^2)*(M^2) - (1/6)*(h^3)*(M^3));
% p3 = zeros(lp,lt);
% tic
% p3(:,1) = pinit;
% for i=2:lt
%     p3(:,i) = A3\(p3(:,i-1));
% end
% t3 = toc;
% P3 = p3';

%%
figure(1)
clf
hold on
h(1) = plot(T, P*ii/params.N,'-ok','DisplayName',['ode45 -> ' num2str(t45,2) 'sec.']);
h(2) = plot(T, P1*ii/params.N,'-xr','DisplayName',['Order 1 -> ' num2str(t1,2) 'sec.']);
h(3) = plot(T, P2*ii/params.N,'-vm','DisplayName',['Order 2 -> ' num2str(t2,2) 'sec.']);
% h(4) = plot(T, P3*ii/params.N,'-*b','DisplayName',['Order 3 -> ' num2str(t3,2) 'sec.']);
hold off
legend(h);
box on
xlabel('Time')
ylabel('Proportion of population infectious')

