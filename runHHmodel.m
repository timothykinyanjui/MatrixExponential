% Set up the data
clearvars; figure

% Load true value calculates by running ODE45 with strict tolerance 1e-13
true = load('TrueValue_Homo','P');
pTrue = true.P; clear true

% Number of people in household
N = 10;

% Set up transmission parameter within the household
b = 0.01;
alpha = 0.663;
beta = b/((N-1)^alpha); %#ok<*PFOUS>

% Transmission between households
% tau = 0.0047;
tau = 0;

% Create the generator matrix
[Q,HHconfig] = SIOld(N); %#ok<*AGROW>

% Generate the data
for i = 1 : length(pTrue(:,1))
    STrue(i,1) = (pTrue(i,:)*HHconfig.dataI(:,1))/N;
    ITrue(i,1) = (pTrue(i,:)*HHconfig.dataI(:,2))/N;
end

% Generate the initial conditions vector
P0 = zeros(length(HHconfig.dataI(:,1)),1); P0(length(HHconfig.dataI(:,1))-1,1) = 1;

% Important shared pameters
% Step size
h = 0.1;
tol = 1e-10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ODE4
% Runge-Kutta order 4 with a fixed time step
f = @(t,x)GenMatrixCalc(Q,beta,tau,HHconfig,x,N)*x;
%h = 10;
timeC = 0:h:365;
% tol = 1e-10;
tic;
P = ode4(f,timeC,P0); %,odeset('RelTol',tol,'NonNegative',1:length(P0)));
timeD(1) = toc;

% Generate the data
for i = 1 : length(P(:,1))
    %E(i,1) = (P(i,:)*HHconfig.dataI(:,2))/N;
    S(i,1) = (P(i,:)*HHconfig.dataI(:,1))/N;
    I(i,1) = (P(i,:)*HHconfig.dataI(:,2))/N;
end
% Check the tolerance at the last time point
tolerance(1) = sqrt(sum((pTrue(end,:)-P(end,:)).^2));
hold on
set(gcf,'WindowStyle','Docked')
subplot(2,3,1)
plot(timeC,S,'b',timeC,I,'r'); hand = legend('S','I'); set(hand,'Box','off')
%xlabel('Time in days')
ylabel('Proportion of individuals')
box off; hold on;
mess = sprintf('ODE4 - %.2f',timeD(1)); title(mess)
% plot(time,S+I,'k--'); set(gca,'YLim',[0 1]) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DA order 1
% DA Method order 1
order  = 1;
%h = 10;
II = eye(length(HHconfig.dataI(:,1)),length(HHconfig.dataI(:,1)));
Mfull = GenMatrixCalc(Q,beta,tau,HHconfig,P0,N);
tic;
pDA_1 = daMethodTime(h,II,Mfull,365,order,P0,beta,tau,HHconfig,N,Q);
timeD(2) = toc;

% Plot
subplot(2,3,2)
hold on
% Generate the data
for i = 1 : length(pDA_1(:,1))
    %E(i,1) = (P(i,:)*HHconfig.dataI(:,2))/N;
    SDA(i,1) = (pDA_1(i,:)*HHconfig.dataI(:,1))/N;
    IDA(i,1) = (pDA_1(i,:)*HHconfig.dataI(:,2))/N;
end
%tolerance(2) = abs(ITrue(end)-IDA(end));
tolerance(2) = sqrt(sum((pTrue(end,:)-pDA_1(end,:)).^2));
hold on
plot(0:h:365,SDA,'b',0:h:365,IDA,'r')
%xlabel('Time in days')
box off; hold on; mess = sprintf('DA order 1 - %.2f',timeD(2)); title(mess)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cheby
% Cheby expansion
%tol = 1e-6;
%h = 10;
timeC = 0:h:365;
tic;
for i = 2:length(timeC)
    pCheb(:,i) = polycheby2(Mfull*timeC(i), P0, tol, 2850, min(diag(Mfull*timeC(i))), max(diag(Mfull*timeC(i))));
    Mfull = GenMatrixCalc(Q,beta,tau,HHconfig,pCheb(:,i),N);
end
timeD(3) = toc;
pCheb(:,1) = P0;
pCheb = pCheb';

% Plot
% Generate the data
for i = 1 : length(pCheb(:,1))
    %E(i,1) = (P(i,:)*HHconfig.dataI(:,2))/N;
    SDC(i,1) = (pCheb(i,:)*HHconfig.dataI(:,1))/N;
    IDC(i,1) = (pCheb(i,:)*HHconfig.dataI(:,2))/N;
end
%tolerance(3) = abs(ITrue(end)-IDC(end));
tolerance(3) = sqrt(sum((pTrue(end,:)-pCheb(end,:)).^2));
subplot(2,3,3); hold on
plot(timeC,SDC,'b',timeC,IDC,'r')
xlabel('Time in days'); ylabel('Proportion of individuals')
box off; mess = sprintf('Cheby Exp - %.2f',timeD(3)); title(mess)
plot(0:h:365,SDC+IDC,'k--'); set(gca,'YLim',[0 1])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Expokit
% Expokit - Krylov Subspace Approximation
%tol = 1e-6;
%h = 10;
timeC = 0:h:365;
Mfull = GenMatrixCalc(Q,beta,tau,HHconfig,P0,N);
tic;
for i = 2:length(timeC)
    pExp(:,i) = mexpv(timeC(i), Mfull, P0, tol);
    Mfull = GenMatrixCalc(Q,beta,tau,HHconfig,pExp(:,i),N);
end
timeD(4) = toc;
pExp(:,1) = P0;
pExp = pExp';
% Plot
% Generate the data
for i = 1 : length(pExp(:,1))
    %E(i,1) = (P(i,:)*HHconfig.dataI(:,2))/N;
    SDEx(i,1) = (pExp(i,:)*HHconfig.dataI(:,1))/N;
    IDEx(i,1) = (pExp(i,:)*HHconfig.dataI(:,2))/N;
end
% tolerance(4) = abs(ITrue(end)-IDEx(end));
tolerance(4) = sqrt(sum((pTrue(end,:)-pExp(end,:)).^2));
subplot(2,3,4); hold on
plot(0:h:365,SDEx,'b',0:h:365,IDEx,'r')
xlabel('Time in days')
box off; mess = sprintf('KSA - %.2f',timeD(4)); title(mess)
plot(0:h:365,SDEx+IDEx,'k--'); set(gca,'YLim',[0 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DA order 1
% DA Method order 2
order  = 2;
%h = 10;
II = eye(length(HHconfig.dataI(:,1)),length(HHconfig.dataI(:,1)));
Mfull = GenMatrixCalc(Q,beta,tau,HHconfig,P0,N);
tic;
pDA_2 = daMethodTime(h,II,Mfull,365,order,P0,beta,tau,HHconfig,N,Q);
pDA_2(pDA_2<0) = 0;
timeD(5) = toc;

% Plot
subplot(2,3,5)
hold on
% Generate the data
for i = 1 : length(pDA_2(:,1))
    %E(i,1) = (P(i,:)*HHconfig.dataI(:,2))/N;
    SDA_2(i,1) = (pDA_2(i,:)*HHconfig.dataI(:,1))/N;
    IDA_2(i,1) = (pDA_2(i,:)*HHconfig.dataI(:,2))/N;
end
%tolerance(2) = abs(ITrue(end)-IDA(end));
tolerance(5) = sqrt(sum((pTrue(end,:)-pDA_2(end,:)).^2));
hold on
plot(0:h:365,SDA_2,'b',0:h:365,IDA_2,'r')
%xlabel('Time in days')
box off; hold on; mess = sprintf('DA order 2 - %.2f',timeD(5)); title(mess)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mohy & Higham et al 2009
% Mohy and Higham new scaling and squaring algorithm for the matrix exponential
%tol = 1e-6;
%h = 10;
timeC = 0:h:365;
Mfull = GenMatrixCalc(Q,beta,tau,HHconfig,P0,N);
tic;
for i = 2:length(timeC)
    pExp_new(:,i) = mexpv_new(timeC(i), Mfull, P0, tol);
    Mfull = GenMatrixCalc(Q,beta,tau,HHconfig,pExp_new(:,i),N);
end
timeD(6) = toc;
pExp_new(:,1) = P0;
pExp_new = pExp_new';
% Plot
% Generate the data
for i = 1 : length(pExp_new(:,1))
    %E(i,1) = (P(i,:)*HHconfig.dataI(:,2))/N;
    SDEx_new(i,1) = (pExp_new(i,:)*HHconfig.dataI(:,1))/N;
    IDEx_new(i,1) = (pExp_new(i,:)*HHconfig.dataI(:,2))/N;
end
% tolerance(4) = abs(ITrue(end)-IDEx(end));
tolerance(6) = sqrt(sum((pTrue(end,:)-pExp_new(end,:)).^2));
subplot(2,3,6); hold on
plot(0:h:365,SDEx_new,'b',0:h:365,IDEx_new,'r')
xlabel('Time in days')
box off; mess = sprintf('Mohy et al - %.2f',timeD(6)); title(mess)
plot(0:h:365,SDEx_new+IDEx_new,'k--'); set(gca,'YLim',[0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot tolerance and accuracy
figure; set(gcf,'WindowStyle','Docked')
gca1 = subplot(1,2,1);
plot(1:6,tolerance,'b*-')
title('Sqrt of Sum of Sqrd differences'); ylabel('Accuracy')
gca2 = subplot(1,2,2);
plot(1:6,timeD,'b*-')
ylabel('Time in seconds')
title('Computational time')
set([gca1, gca2],'XtickLabel',{'ODE4','DA-1','Cheby','KSA','DA-2','Mohy et al'},'XLim',[1 6],'Box','off')
