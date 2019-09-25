function P = daMethodTime(h,II,M,T,order,PP0,beta,tau,gamma,HHconfig,N,Q)

% Function solves using the Backward Euler method
% Written by Tim on 7th Aug 2015
% University of Manchester

% Set time vector
time = 0:h:T;

% Solve using Euler backward of order order
P(:,1) = PP0;

if order == 1
    
    for i = 2:length(time)
        M = GenMatrixCalc(Q,beta,tau,gamma,HHconfig,P(:,i-1),N);
        Msum = II - h*M;
        P(:,i) = Msum\P(:,i-1);
    end
    
elseif order == 2
    
    for i = 2:length(time)
        M = GenMatrixCalc(Q,beta,tau,gamma,HHconfig,P(:,i-1),N);
        Msum = II - h*M + ((h^2)/prod(1:2))*(M*M);
        P(:,i) = Msum\P(:,i-1);
    end
    
elseif order == 3
    
    for i = 2:length(time)
        M = GenMatrixCalc(Q,beta,tau,gamma,HHconfig,P(:,i-1),N);
        Msum = II - h*M + ((h^2)/prod(1:2))*(M*M) - ((h^3)/prod(1:3))*(M*M*M);
        P(:,i) = Msum\P(:,i-1);
    end
    
else
    
    error('Can only accomodate order 1,2 and 3')
    
end

P = P';

return