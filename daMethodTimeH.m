function P = daMethodTimeH(h,II,T,PP0,beta,tau,HHconfig,N,Q)

% Function solves using the Backward Euler method
% Written by Tim on 7th Aug 2015
% University of Manchester

% Set time vector
time = 0:h:T;

% Solve using Euler backward of order order
P(:,1) = PP0;

% Loop to solve the system
for i = 2:length(time)
    
    % Calculate M
    M = GenMatrixCalc(Q,beta,tau,HHconfig,P(:,i-1),N);
    
    % Derivative of M
    Mdot = tau*((HHconfig.dataI(:,2)'*(M*P(:,i-1)))/N)*Q.QC;
    
    % Second derivative of M
    dMP = Mdot*P(:,i-1) + (M^2)*P(:,i-1);
    Mdotdot = tau*((HHconfig.dataI(:,2)'*dMP)/N)*Q.QC;
    
    % Do the Taylor expansion
    %Msum = II + h*M + (h^2)*((M^2)+Mdot)/2 + (h^3)*(Mdotdot + 2*(Mdot*M) + (M*Mdot) + M^3)/6;
    Msum = II + h*M + (h^2)*((M^2)+Mdot)/2 + (h^3)*((M^3) + 6*(M*Mdot)/2 + 6*(Mdotdot + 2*(Mdot*M) - 2*(M*Mdot)))/6; % Product form
    
    % The forward step
    P(:,i) = Msum*P(:,i-1);
    
end

P = P';

% End
return