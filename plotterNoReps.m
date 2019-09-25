clear
% Load the data
%load refined_N_270516.mat
load testRun
D = D(1:4);

% Tolerance
tolerance = log(tolerance);
meanTol = tolerance;
%D = matsize;

% Do the plotting
figure; set(gcf,'WindowStyle','docked')
subplot(1,2,1)
% [hl, hp] = boundedline(D,meanTol(:,1),[lb(:,1) ub(:,1)],'r','alpha',D,meanTol(:,2),[lb(:,2) ub(:,2)],'b','alpha',D,meanTol(:,3),[lb(:,3) ub(:,3)],'g','alpha',D,meanTol(:,4),[lb(:,4) ub(:,4)],'k','alpha',D,meanTol(:,5),[lb(:,5) ub(:,5)],'m','alpha',D,meanTol(:,6),[lb(:,6) ub(:,6)],'c','alpha',D,meanTol(:,7),[lb(:,7) ub(:,7)],'y','alpha');
% outlinebounds(hl,hp);
plot(D,meanTol(:,1),'r',D,meanTol(:,2),'b',D,meanTol(:,3),'g',D,meanTol(:,4),'k',D,meanTol(:,5),'m',D,meanTol(:,6),'c',D,meanTol(:,7),'y',D,meanTol(:,8),'b--');
title('SSSD'); ylabel('Accuracy (SSSD)'); xlabel('System size')
hand = legend('ODE4','DA-1','Cheby','KSA','DA-2','Mohy','DA-3','Mohy2010','Location','Best'); set(hand,'Box','off')
box off

% Time
timeD = log(timeD);
meanTime = timeD;

% Do the plotting
subplot(1,2,2)
plot(D,meanTime(:,1),'r',D,meanTime(:,2),'b',D,meanTime(:,3),'g',D,meanTime(:,4),'k',D,meanTime(:,5),'m',D,meanTime(:,6),'c',D,meanTime(:,7),'y',D,meanTime(:,8),'b--');
title('Computational time'); ylabel('Time in log seconds'); xlabel('System size')
box off