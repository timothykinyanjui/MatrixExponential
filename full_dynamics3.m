% Script to integrate the full system - plots simple infected curve rather
% than complex figure in main paper

clear

params.beta = 3;
params.N = 200;
params.g = 1;
params.mu = 0;
params.alpha = 0;
tmax = 3;

tol = 0.01; % TOL parameter for variable timestepper 
tol2 = 0.01;
tol3 = 0.01;
h = 0.01;

% [M,ss,ii] = SIRS(params.N,params.beta/(params.N-1),params.g,params.mu,0,[0 params.N],[0 params.N]);
% lp = length(ii);
% params.ii = ii;
% pinit = zeros(length(ii),1);
% pinit(ii==1 & ss==(params.N-1))=1;
%BELOW USES REACTION COUNTS --> DA TRIANGULAR REPRESENTATION
m = SIR_DA(params.N-1,1);
M = input_params_tran(m,params.beta,params.g);
lp = length(M);
pinit = zeros(lp,1);
pinit(1) = 1;
ii = 0:params.N;
TM = create_BImat(m);


MNN = M - diag(diag(M));
Msq = M^2;
MsqNN = Msq - diag(diag(Msq));
Mp3 = M^3;
Mp3NN = Mp3 - diag(diag(Mp3));
MODE = max(max(1/2*h*MsqNN - MNN - 1/6*h^2*Mp3));
%ho = min(min(2*MNN./MsqNN))

Qdiff = @(t,p) M*p;

tran = 0:h:tmax;
lt = length(tran);

tic
[T, P] = ode45(Qdiff, tran, pinit);
t45 = toc;

opts.LT = true;

A1 = speye(lp) - h*M;
p1 = zeros(lp,lt);
tic
p1(:,1) = pinit;
n1 = 0;
for i=2:lt
    p1(:,i) = A1\(p1(:,i-1));
    P1 = p1';
    n1(i) = norm(P1(i,:)-P(i,:));
end
t1 = toc;
P1 = p1';

A2 = speye(lp) - h*M + (1/2)*(h^2)*(M^2);
p2 = zeros(lp,lt);
tic
p2(:,1) = pinit;
n2 = 0;
for i=2:lt
    p2(:,i) = A2\(p2(:,i-1));
    if min(p2(:,i)) < 0
        p2(find(p2(:,i)<0),i) = zeros(1,length(find(p2(:,i)<0)));
        p2(:,i) = p2(:,i)/sum(p2(:,i));
    end
    P2 = p2';
    n2(i) = norm(P2(i,:)-P(i,:));
end
t2 = toc;
P2 = p2';

A3 = speye(lp) - h*M + ...
    (1/2)*(h^2)*(M^2) - (1/6)*(h^3)*(M^3);
p3 = zeros(lp,lt);
tic
p3(:,1) = pinit;
n3 = 0;
for i=2:lt
    p3(:,i) = A3\(p3(:,i-1));
    if min(p3(:,i)) < 0
        p3(find(p3(:,i)<0),i) = zeros(1,length(find(p3(:,i)<0)));
        p3(:,i) = p3(:,i)/sum(p3(:,i));
    end
    P3 = p3';
    n3(i) = norm(P3(i,:)-P(i,:));
end
t3 = toc;
P3 = p3';

%Richardson Extrapolation with Fixed Step Size
A1h = speye(lp) - h/2*M;
p1h = zeros(lp,lt);
tic
p1h(:,1) = pinit;
nRF = 0;
for i=2:lt
    p1h1(:,i) = A1\(p1h(:,i-1));
    p1ht(:,i) = A1h\(p1h(:,i-1));
    p1h2(:,i) = A1h\(p1ht(:,i));
    p1h(:,i) = 2*p1h2(:,i) - p1h1(:,i);
    if min(p1h(:,i)) < 0
        p1h(:,i) = p1h2(:,i);
    end
    P1h = p1h';
    nRF(i) = norm(P1h(i,:)-P(i,:));
end
t1h = toc;
P1h = p1h';

%Richardson Extrapolation with Variable Step Size
tnow = 0;
tot_error = zeros(lp,1);
h = tol;
p1RV = zeros(lp,2);
p1RV1 = zeros(lp,2);
p1RV2 = zeros(lp,2);
tic
p1RV(:,1) = pinit;
p1RV1(:,1) = pinit;
p1RV2(:,1) = pinit;
i = 2;
while tnow(end) < tmax
    A1 = speye(lp) - h*M;
    A1h = speye(lp) - h/2*M;
    p1RV1(:,i) = A1\(p1RV(:,i-1));
    p1RV2t(:,i) = A1h\(p1RV(:,i-1));
    p1RV2(:,i) = A1h\(p1RV2t(:,i));
    error_loc = (p1RV2(:,i)-p1RV1(:,i));
    normError = norm(error_loc,1);
    tol_loc = tol; 
    if normError <= tol_loc %accept step, increase step size
        tot_error = tot_error + error_loc;
        % set new probability via Richardson Extrapolation
        %
        p1RV(:,i) = 2*p1RV2(:,i) - p1RV1(:,i); 
        if min(p1RV(:,i)) < 0
            p1RV(:,i) = p1RV2(:,i);
        end
        % update time
        tnow(i)  = tnow(i-1) + h; 
        % adaptively choose stepsize
        factor = min(10 , .9*sqrt(tol_loc/normError)); 
        h = h*factor;
        % round stepsize for numerical stability
        s = 10^(floor(log10(h))-1);
        h = min(ceil(h/s) * s, tmax - tnow(end));
        i = i+1;
     else % do not accept step, decrease step size       
        factor = max(1/5 , .9*sqrt(tol_loc/normError));
        h = h*factor;% set stepsize
        % round stepsize for numerical stability
        s = 10^(floor(log10(h))-1);
        h = min(ceil(h/s) * s, tmax - tnow(end));       
    end
end
t1RV = toc;
P1RV = p1RV';
[TRV, PRVc] = ode45(Qdiff, tnow, pinit);
for i = 1:length(tnow)
    nRV(i) = norm(P1RV(i,:) - PRVc(i,:));
end



%Our Order 2 with Variable Step Size
tnow2 = 0;
tot_error = zeros(lp,1);
h = tol2;
p1OV = zeros(lp,2);
p1OV1 = zeros(lp,2);
p1OV2 = zeros(lp,2);
tic
p1OV(:,1) = pinit;
p1OV1(:,1) = pinit;
p1OV2(:,1) = pinit;
i = 2;
Ms = M^2;
while tnow2(end) < tmax
    %A1 = speye(lp) - h*M;
    A2 = speye(lp) - h*M + (1/2)*(h^2)*(Ms);
    %A1h = speye(lp) - h/2*M;
    A2h = speye(lp) - h/2*M + (1/2)*((h/2)^2)*(Ms);
    p1OV1(:,i) = A2\(p1OV(:,i-1));
    p1OV2t(:,i) = A2h\(p1OV(:,i-1));
    p1OV2(:,i) = A2h\(p1OV2t(:,i));
    error_loc = (p1OV2(:,i)-p1OV1(:,i));
    normError = norm(error_loc,1);
    tol_loc = tol2; 
    if normError <= tol_loc %accept step, increase step size
        tot_error = tot_error + error_loc;
        % set new probability
        %
        p1OV(:,i) = p1OV2(:,i); 
        if min(p1OV(:,i)) < 0
            p1OV(find(p1OV(:,i)<0),i) = zeros(1,length(find(p1OV(:,i)<0)));
            p1OV(:,i) = p1OV(:,i)/sum(p1OV(:,i));
        end
        % update time
        tnow2(i)  = tnow2(i-1) + h; 
        % adaptively choose stepsize
        factor = min(10,.9*sqrt(tol_loc/normError)); 
        h = h*factor;
        % round stepsize for numerical stability
        s = 10^(floor(log10(h))-1);
        h = min(ceil(h/s) * s, tmax - tnow2(end));
        i = i+1;
     else % do not accept step, decrease step size       
        factor = max(1/5,.9*sqrt(tol_loc/normError));
        h = h*factor;% set stepsize
        % round stepsize for numerical stability
        s = 10^(floor(log10(h))-1);
        h = min(ceil(h/s) * s, tmax - tnow2(end));       
    end
end
t1OV = toc;
P1OV = p1OV';
[TRV, POVc] = ode45(Qdiff, tnow2, pinit);
for i = 1:length(tnow2)
    nOV2(i) = norm(P1OV(i,:) - POVc(i,:));
end



%Our Order 3 with Variable Step Size
tnow3 = 0;
tot_error = zeros(lp,1);
h = tol3;
p2OV = zeros(lp,2);
p1OV1 = zeros(lp,2);
p1OV2 = zeros(lp,2);
tic
p2OV(:,1) = pinit;
p1OV1(:,1) = pinit;
p1OV2t(:,1) = pinit;
p1OV2(:,1) = pinit;
i = 2;
Mp3 = M^3;
while tnow3(end) < tmax
    %A1 = speye(lp) - h*M;
    A3 = speye(lp) - h*M + ...
    (1/2)*(h^2)*(Ms) - (1/6)*(h^3)*(Mp3);
    A3h = speye(lp) - h/2*M + ...
    (1/2)*((h/2)^2)*(Ms) - (1/6)*((h/2)^3)*(Mp3);
    p1OV1(:,i) = A3\(p2OV(:,i-1));
    p1OV2t(:,i) = A3h\(p2OV(:,i-1));
    p1OV2(:,i) = A3h\(p1OV2t(:,i));
    error_loc = (p1OV2(:,i)-p1OV1(:,i));
    normError = norm(error_loc,1);
    tol_loc = tol3; 
    if normError <= tol_loc %accept step, increase step size
        tot_error = tot_error + error_loc;
        % set new probability
        %
        p2OV(:,i) = p1OV2(:,i); 
        if min(p2OV(:,i)) < 0
            p2OV(find(p2OV(:,i)<0),i) = zeros(1,length(find(p2OV(:,i)<0)));
            p2OV(:,i) = p2OV(:,i)/sum(p2OV(:,i));
        end
        % update time
        tnow3(i)  = tnow3(i-1) + h; 
        % adaptively choose stepsize
        factor = min(10,.9*sqrt(tol_loc/normError)); 
        h = h*factor;
        % round stepsize for numerical stability
        s = 10^(floor(log10(h))-1);
        h = min(ceil(h/s) * s, tmax - tnow3(end));
        i = i+1;
     else % do not accept step, decrease step size       
        factor = max(1/5,.9*sqrt(tol_loc/normError));
        h = h*factor;% set stepsize
        % round stepsize for numerical stability
        s = 10^(floor(log10(h))-1);
        h = min(ceil(h/s) * s, tmax - tnow3(end));       
    end
end
t2OV = toc;
P2OV = p2OV';
[TRV, POVc] = ode45(Qdiff, tnow3, pinit);
for i = 1:length(tnow3)
    nOV3(i) = norm(P2OV(i,:) - POVc(i,:));
end



figure(1)
clf
hold on
% h(1) = plot(T, P*ii/params.N,'-ok','DisplayName',['ode45 -> ' num2str(t45,2) 'sec.']);
% h(2) = plot(T, P1*ii/params.N,'-xr','DisplayName',['Order 1 -> ' num2str(t1,2) 'sec.']);
% h(3) = plot(T, P2*ii/params.N,'-vm','DisplayName',['Order 2 -> ' num2str(t2,2) 'sec.']);
% h(4) = plot(T, P3*ii/params.N,'-*b','DisplayName',['Order 3 -> ' num2str(t3,2) 'sec.']);
h(1) = plot(T, n1,'-xr','DisplayName',['Order 1']);
h(2) = plot(T, n2,'-vm','DisplayName',['Order 2']);
h(3) = plot(T, n3,'-*b','DisplayName',['Order 3']);
h(4) = plot(T, nRF,'-.k','DisplayName',['Richardson']);
h(5) = plot(tnow, nRV,'-xg','DisplayName',['Richardson variable']);
h(6) = plot(tnow2, nOV2,'-xy','DisplayName',['Order 2 variable']);
h(7) = plot(tnow3, nOV3,'-or','DisplayName',['Order 3 variable']);
hold off
legend(h);
box on
xlabel('Time')
ylabel('Proportion of population infectious')
