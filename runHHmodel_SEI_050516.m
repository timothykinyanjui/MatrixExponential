% Set up the data
clearvars; figure

% Load true value calculates by running ODE45 with strict tolerance 1e-13
true = load('TrueValue_Homo','P');
pTrue = true.P; clear true

% Number of people in household - Changing N will change the system size
N = 10;

% Set up transmission parameter within the household
b = 0.1;
alpha = 0.663;
beta = b/((N-1)^alpha); %#ok<*PFOUS>
gamma = 0.025;

% Transmission between households
% tau = 0.0047;
tau = 0;

% Create the generator matrices
[Q,HHconfig] = SEI(N); %#ok<*AGROW>

% Generate the data
for i = 1 : length(pTrue(:,1))
    STrue(i,1) = (pTrue(i,:)*HHconfig.dataI(:,1))/N;
    ITrue(i,1) = (pTrue(i,:)*HHconfig.dataI(:,3))/N;
    ETrue(i,1) = (pTrue(i,:)*HHconfig.dataI(:,2))/N;
end

% Generate the initial conditions vector
tempI = find(HHconfig.dataI(:,3)==1); tempS = find(HHconfig.dataI(:,1)==N-1);
pos = intersect(tempI,tempS);
P0 = zeros(length(HHconfig.dataI(:,1)),1); P0(pos,1) = 1;

% Important shared pameters
% Step size
h = 5;
tol = 1e-8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ODE4
% Runge-Kutta order 4 with a fixed time step
f = @(t,x)GenMatrixCalc(Q,beta,tau,gamma,HHconfig,x,N)*x;
%h = 10;
timeC = 0:h:365;
% tol = 1e-10;
tic;
P = ode4(f,timeC,P0); %,odeset('RelTol',tol,'NonNegative',1:length(P0)));
timeD(1) = toc;

% Generate the data
for i = 1 : length(P(:,1))
    E(i,1) = (P(i,:)*HHconfig.dataI(:,2))/N;
    S(i,1) = (P(i,:)*HHconfig.dataI(:,1))/N;
    I(i,1) = (P(i,:)*HHconfig.dataI(:,3))/N;
end
% Check the tolerance at the last time point
%tolerance(1) = sqrt(sum((pTrue(end,:)-P(end,:)).^2));
pODE4 = abs(P(end,:)); % Renormalise this
pODE4 = pODE4/sum(pODE4);
tolerance(1) = kl_Div(pTrue(end,:),pODE4);
hold on
set(gcf,'WindowStyle','Docked')
subplot(2,4,1)
plot(timeC,S,'b',timeC,I,'r',timeC,E,'k'); hand = legend('S','I','E'); set(hand,'Box','off')
%xlabel('Time in days')
ylabel('Proportion')
box off; hold on;
mess = sprintf('ODE4 - %.2f',timeD(1)); title(mess)
plot(timeC,S+I+E,'k--'); set(gca,'YLim',[0 1]) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DA order 1
% DA Method order 1
order  = 1;
%h = 10;
II = eye(length(HHconfig.dataI(:,1)),length(HHconfig.dataI(:,1)));
Mfull = GenMatrixCalc(Q,beta,tau,gamma,HHconfig,P0,N);
tic;
pDA_1 = daMethodTime(h,II,Mfull,365,order,P0,beta,tau,gamma,HHconfig,N,Q);
timeD(2) = toc;

% Plot
subplot(2,4,2)
hold on
% Generate the data
for i = 1 : length(pDA_1(:,1))
    EDA(i,1) = (pDA_1(i,:)*HHconfig.dataI(:,2))/N;
    SDA(i,1) = (pDA_1(i,:)*HHconfig.dataI(:,1))/N;
    IDA(i,1) = (pDA_1(i,:)*HHconfig.dataI(:,3))/N;
end
%tolerance(2) = sqrt(sum((pTrue(end,:)-pDA_1(end,:)).^2));
tolerance(2) = kl_Div(pTrue(end,:),pDA_1(end,:));
hold on
plot(0:h:365,SDA,'b',0:h:365,IDA,'r',0:h:365,EDA,'k')
%xlabel('Time in days')
box off; hold on; mess = sprintf('DA order 1 - %.2f',timeD(2)); title(mess)
plot(0:h:365,SDA+IDA+EDA,'k--'); set(gca,'YLim',[0 1]) 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cheby
% Cheby expansion
%tol = 1e-6;
%h = 10;
timeC = 0:h:365;

tic;
for i = 2:length(timeC)
    Mfullall{i} = Mfull;
    pCheb(:,i) = polycheby2(Mfull*timeC(i), P0, tol, 2850, min(diag(Mfull*timeC(i))), max(diag(Mfull*timeC(i))));
    Mfull = GenMatrixCalc(Q,beta,tau,gamma,HHconfig,pCheb(:,i),N); 
end
timeD(3) = toc;

% This part corrects to ensure result is a probability vector
%pCheb(pCheb<0) = 0;
pCheb(pCheb>1) = 1;
pCheb = abs(pCheb);
pCheb(:,1) = P0;
pCheb = pCheb';
for i = 1:length(pCheb(:,1))
    pCheb(i,:) = pCheb(i,:)/sum(pCheb(i,:));
end

% Plot
% Generate the data
for i = 1 : length(pCheb(:,1))
    EDC(i,1) = (pCheb(i,:)*HHconfig.dataI(:,2))/N;
    SDC(i,1) = (pCheb(i,:)*HHconfig.dataI(:,1))/N;
    IDC(i,1) = (pCheb(i,:)*HHconfig.dataI(:,3))/N;
end
% tolerance(3) = sqrt(sum((pTrue(end,:)-pCheb(end,:)).^2));
tolerance(3) = kl_Div(pTrue(end,:),pCheb(end,:));
subplot(2,4,3); hold on
plot(timeC,SDC,'b',timeC,IDC,'r',timeC,EDC,'k')
box off; mess = sprintf('Cheby Exp - %.2f',timeD(3)); title(mess)
plot(0:h:365,SDC+IDC+EDC,'k--'); set(gca,'YLim',[0 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Expokit
% Expokit - Krylov Subspace Approximation
%tol = 1e-6;
%h = 10;
timeC = 0:h:365;
Mfull = GenMatrixCalc(Q,beta,tau,gamma,HHconfig,P0,N);
tic;
for i = 2:length(timeC)
    pExp(:,i) = mexpv(timeC(i), Mfull, P0, tol);
    Mfull = GenMatrixCalc(Q,beta,tau,gamma,HHconfig,pExp(:,i),N);
end
timeD(4) = toc;
pExp(:,1) = P0;
pExp = pExp';
% Plot
% Generate the data
for i = 1 : length(pExp(:,1))
    EDEx(i,1) = (pExp(i,:)*HHconfig.dataI(:,2))/N;
    SDEx(i,1) = (pExp(i,:)*HHconfig.dataI(:,1))/N;
    IDEx(i,1) = (pExp(i,:)*HHconfig.dataI(:,3))/N;
end
% tolerance(4) = sqrt(sum((pTrue(end,:)-pExp(end,:)).^2));
tolerance(4) = kl_Div(pTrue(end,:),pExp(end,:));
subplot(2,4,4); hold on
plot(0:h:365,SDEx,'b',0:h:365,IDEx,'r',0:h:365,EDEx,'k')
%xlabel('Time in days')
box off; mess = sprintf('KSA - %.2f',timeD(4)); title(mess)
plot(0:h:365,SDEx+IDEx+EDEx,'k--'); set(gca,'YLim',[0 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DA order 2
% DA Method order 2
order  = 2;
%h = 10;
II = eye(length(HHconfig.dataI(:,1)),length(HHconfig.dataI(:,1)));
Mfull = GenMatrixCalc(Q,beta,tau,gamma,HHconfig,P0,N);
tic;
pDA_2 = daMethodTime(h,II,Mfull,365,order,P0,beta,tau,gamma,HHconfig,N,Q);
pDA_2(pDA_2<0) = 0;
timeD(5) = toc;

% Plot
subplot(2,4,5)
hold on
% Generate the data
for i = 1 : length(pDA_2(:,1))
    EDA_2(i,1) = (pDA_2(i,:)*HHconfig.dataI(:,2))/N;
    SDA_2(i,1) = (pDA_2(i,:)*HHconfig.dataI(:,1))/N;
    IDA_2(i,1) = (pDA_2(i,:)*HHconfig.dataI(:,3))/N;
end
% tolerance(5) = sqrt(sum((pTrue(end,:)-pDA_2(end,:)).^2));
tolerance(5) = kl_Div(pTrue(end,:),pDA_2(end,:));
hold on
plot(0:h:365,SDA_2,'b',0:h:365,IDA_2,'r',0:h:365,EDA_2,'k')
xlabel('Time in days'); ylabel('Proportion')
box off; hold on; mess = sprintf('DA order 2 - %.2f',timeD(5)); title(mess)
plot(0:h:365,SDA_2+IDA_2+EDA_2,'k--'); set(gca,'YLim',[0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mohy & Higham et al 2009
% Mohy and Higham new scaling and squaring algorithm for the matrix exponential
%tol = 1e-6;
%h = 10;
timeC = 0:h:365;
Mfull = GenMatrixCalc(Q,beta,tau,gamma,HHconfig,P0,N);
tic;
for i = 2:length(timeC)
    pExp_new(:,i) = mexpv_new(timeC(i), Mfull, P0, tol);
    Mfull = GenMatrixCalc(Q,beta,tau,gamma,HHconfig,pExp_new(:,i),N);
end
timeD(6) = toc;
pExp_new(:,1) = P0;
pExp_new = pExp_new';
% Plot
% Generate the data
for i = 1 : length(pExp_new(:,1))
    EDEx_new(i,1) = (pExp_new(i,:)*HHconfig.dataI(:,2))/N;
    SDEx_new(i,1) = (pExp_new(i,:)*HHconfig.dataI(:,1))/N;
    IDEx_new(i,1) = (pExp_new(i,:)*HHconfig.dataI(:,3))/N;
end
% tolerance(6) = sqrt(sum((pTrue(end,:)-pExp_new(end,:)).^2));
tolerance(6) = kl_Div(pTrue(end,:),pExp_new(end,:));
subplot(2,4,6); hold on
plot(0:h:365,SDEx_new,'b',0:h:365,IDEx_new,'r',0:h:365,EDEx_new,'k')
xlabel('Time in days')
box off; mess = sprintf('Mohy et al - %.2f',timeD(6)); title(mess)
plot(0:h:365,SDEx_new+IDEx_new+EDEx_new,'k--'); set(gca,'YLim',[0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DA order 3
% This is the backward Euler order 3
order  = 3;
%h = 10;
II = eye(length(HHconfig.dataI(:,1)),length(HHconfig.dataI(:,1)));
Mfull = GenMatrixCalc(Q,beta,tau,gamma,HHconfig,P0,N);
tic;
pDA_3 = daMethodTime(h,II,Mfull,365,order,P0,beta,tau,gamma,HHconfig,N,Q);
timeD(7) = toc;

% Plot
subplot(2,4,7)
hold on
% Generate the data
for i = 1 : length(pDA_3(:,1))
    EDA_3(i,1) = (pDA_3(i,:)*HHconfig.dataI(:,2))/N;
    SDA_3(i,1) = (pDA_3(i,:)*HHconfig.dataI(:,1))/N;
    IDA_3(i,1) = (pDA_3(i,:)*HHconfig.dataI(:,3))/N;
end
% tolerance(7) = sqrt(sum((pTrue(end,:)-pDA_3(end,:)).^2));
tolerance(7) = kl_Div(pTrue(end,:),pDA_3(end,:));
hold on
plot(0:h:365,SDA_3,'b',0:h:365,IDA_3,'r',0:h:365,EDA_3,'k')
xlabel('Time in days')
box off; hold on; mess = sprintf('DA order 3 - %.2f',timeD(7)); title(mess)
plot(0:h:365,SDA_3+IDA_3+EDA_3,'k--'); set(gca,'YLim',[0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mohy & Higham 2010
timeC = 0:h:365;
Mfull = GenMatrixCalc(Q,beta,tau,gamma,HHconfig,P0,N); % Calculates the Q matrix
tic;
for i = 2:length(timeC)
    [pNew(:,i),s,m,mv,mvd,unA] = expmv_2010(timeC(i),Mfull,P0,[],'single');     %#ok<*SAGROW>
end
timeD(8) = toc;
pNew(:,1) = P0;
pNew = pNew';
% This part corrects to make probability vector
for i = 1:length(pNew(:,1))
    pNew(i,:) = pNew(i,:)/sum(pNew(i,:));
end

% Plot
for i = 1 : length(pNew(:,1))
    Enew(i,1) = (pNew(i,:)*HHconfig.dataI(:,2))/N;
    Snew(i,1) = (pNew(i,:)*HHconfig.dataI(:,1))/N;
    Inew(i,1) = (pNew(i,:)*HHconfig.dataI(:,3))/N;
end
subplot(2,4,8); hold on
plot(timeC,Snew,'b',timeC,Enew,'k',timeC,Inew,'r')
box off; mess = sprintf('Mohy & Higham (2010) - %.2f',timeD(8)); title(mess)
plot(0:h:365,Snew+Inew+Enew,'k--'); set(gca,'YLim',[0 1])
xlabel('Time in days')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot tolerance and accuracy

figure; set(gcf,'WindowStyle','Docked')
gca1 = subplot(1,2,1);
plot(1:7,log(tolerance),'b*-')
title('K-L Divergence'); ylabel('log Accuracy/Error')
gca2 = subplot(1,2,2);
plot(1:8,log(timeD),'b*-')
ylabel('Time in log seconds')
title('Computational time')
set([gca1, gca2],'XtickLabel',{'ODE4','DA-1','Cheby','KSA','DA-2','Mohy et al','DA-3'},'XLim',[1 7],'Xtick',1:7,'Box','off','XTickLabelRotation',45)
