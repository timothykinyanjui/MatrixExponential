clearvars
% Load the data
load ModelRunsReplicates_h_SSSD


% Tolerance
tolerance = log(tolerance);
meanTol = mean(tolerance,3);

% Calculate the prctile
Y = prctile(tolerance,[2.5 97.5],3); 
lb = abs(meanTol - Y(:,:,1));
ub = abs(meanTol - Y(:,:,2));

% Do the plotting
figure; set(gcf,'WindowStyle','docked')
subplot(1,2,1)
% [hl, hp] = boundedline(D,meanTol(:,1),[lb(:,1) ub(:,1)],'r','alpha',D,meanTol(:,2),[lb(:,2) ub(:,2)],'b','alpha',D,meanTol(:,3),[lb(:,3) ub(:,3)],'g','alpha',D,meanTol(:,4),[lb(:,4) ub(:,4)],'k','alpha',D,meanTol(:,5),[lb(:,5) ub(:,5)],'m','alpha',D,meanTol(:,6),[lb(:,6) ub(:,6)],'c','alpha',D,meanTol(:,7),[lb(:,7) ub(:,7)],'y','alpha');
% outlinebounds(hl,hp);
hand = plot(D,meanTol(:,1),'r',D,meanTol(:,2),'b',D,meanTol(:,3),'g',D,meanTol(:,4),'k',D,meanTol(:,5),'m',D,meanTol(:,6),'c',D,meanTol(:,7),'y',D,meanTol(:,8),'b--');
title('Accuracy - SSSD'); ylabel('log(SSSD)'); xlabel('Step size in days')
set(hand,'LineWidth',1.5)
hand = legend('ODE4','DA-1','Cheby','KSA','DA-2','Mohy','DA-3','Mohy2010','Location','Best'); set(hand,'Box','off')
set(gca,'XTick',0:50:400)
box off; hand = text(350,1.5,'A'); set(hand,'FontWeight','bold')

% Time
timeD = log(timeD);
meanTime = mean(timeD,3);

% Calculate the prctile
Y = prctile(timeD,[2.5 97.5],3); 
lb = abs(meanTime - Y(:,:,1));
ub = abs(meanTime - Y(:,:,2));

subplot(1,2,2)
[hl, hp] = boundedline(D,meanTime(:,1),[lb(:,1) ub(:,1)],'r','alpha',D,meanTime(:,2),[lb(:,2) ub(:,2)],'b','alpha',D,meanTime(:,3),[lb(:,3) ub(:,3)],'g','alpha',D,meanTime(:,4),[lb(:,4) ub(:,4)],'k','alpha',D,meanTime(:,5),[lb(:,5) ub(:,5)],'m','alpha',D,meanTime(:,6),[lb(:,6) ub(:,6)],'c','alpha',D,meanTime(:,7),[lb(:,7) ub(:,7)],'y',D,meanTime(:,8),[lb(:,8) ub(:,8)],'b--','alpha');
%outlinebounds(hl,hp);
title('Computational time'); ylabel('Time (log seconds)'); xlabel('Step size in days')
box off; set(gca,'XTick',0:50:400); hand = text(350,3.5,'B'); set(hand,'FontWeight','bold')