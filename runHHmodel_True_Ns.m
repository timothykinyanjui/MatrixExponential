clearvars

% Set up transmission parameter within the household
b = 0.1;
alpha = 0.663;
gamma = 0.025;
h = 1;


% Transmission between households
% tau = 0.0047;
tau = 0;
% NN = [2:2:60 70 80 90];
NN = [30 99]; % Household size

parfor ii = 1:length(NN)
    
    N = NN(ii);
    beta = b/((N-1)^alpha);
    
    % Create the generator matrix
    [Q,HHconfig] = SEI(N); %#ok<*AGROW>
    
    % Generate the initial conditions vector
    % Generate the initial conditions vector
    tempI = find(HHconfig.dataI(:,3)==1); tempS = find(HHconfig.dataI(:,1)==N-1);
    pos = intersect(tempI,tempS);
    P0 = zeros(length(HHconfig.dataI(:,1)),1); P0(pos,1) = 1;
    
    % Tolerance
    tol = 1e-12; % Similar to what I have for many h's
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ODE4
    % Runge-Kutta order 4 with a fixed time step
    f = @(t,x)GenMatrixCalc(Q,beta,tau,gamma,HHconfig,x,N)*x;
    timeC = 0:h:360;
    tic;
    [tt, PP] = ode45(f,timeC,P0,odeset('RelTol',tol,'NonNegative',1:length(P0)));%,'NonNegative',1:length(P0)));
    P(ii).PP = PP(end,:);
    
end

% Save
save TrueValue_Homo_Ns