function Mfull = GenMatrixCalc(Q,beta,tau,gamma,HHconfig,P0,N)

% Function calculates the generator matrix for the time inhomogeneous case
% for the SI model.
% Written by Tim on 16 Nov 2015
% At University of Adelaide, SA

% Set up relevant stuff
Q.inf = beta*Q.inf;

M = Q.inf + tau*((HHconfig.dataI(:,2)'*P0)/N)*Q.QC + gamma*Q.prog;

Mfull = M'; % Transpose the matrix

% End
return