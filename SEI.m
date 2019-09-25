function [Q, HHconfig] = SEI(N)

% Function written to generate the transition rates matrix for household
% model with SI disease structure
%
% Input:
%   N: Household size which we hold constant in the simulations
%
% Output:
%   Q: Transition matrix
%   dataI: Contains the household configurations
%
% Function written by Tim Kinyanjui on 16th Sept 2015
% University of Manchester


% List the household configurations exhaustively i.e. all the possible
% household configurations and exclude the ones that do not retain the
% household size


% Dimension of the matrix. This only works for 3 epidemiological classes. I
% don't know what will happen when they are not 3, but will find out
% later.
dim = sum(1:(N+1));

Q.inf = zeros(dim,dim);
Q.QC = zeros(dim,dim);
Q.prog = zeros(dim,dim);

dataI = zeros(dim,3);

m = 0;

% Transforms two integers to locations where we will store the variables
for ss=0:N
    
    for ii=0:(N-ss)
        
        m = m+1;
        
        I{ss+1,ii+1} = m;  %#ok<AGROW>
        
    end
    
end

% Counter for susceptibles
for ss=0:N
    
    for ii=0:(N-ss)
        
        % If susceptibles are less than N then there are infecteds and therefore
        % infection within the household can happen
        if  ss >= 1 && (N-ss-ii) >=1
            
            % Rate of within household infection
            % Q.inf(I{ss+1,ii+1},I{ss}) = betaHH*ss*(N-ss)/(N-1);
            Q.inf(I{ss+1,ii+1},I{ss+1-1,ii+1+1}) = ss*(N-ss-ii);
            
        end
        
        % Infection from outside of the household
        % Infection will only occur if number of susceptibles is greater than 0
        if ss >= 1
            
            % Infection from outside of the household
            % Q.QC(I{ss+1},I{ss+1-1}) = ss;
            Q.QC(I{ss+1,ii+1},I{ss+1-1,ii+1+1}) = ss;
            
        end
        
        % For progression to infection to occur
        if ii >= 1
            
            % Progression to infection
            Q.prog(I{ss+1,ii+1},I{ss+1,ii+1-1}) = ii;
            
        end
        
        % Store the relevant indices to help identify the household
        % configurations outside of this function.
        
        dataI(I{ss+1,ii+1},:) = [ss, ii, N-ss-ii];
        
    end
    
end

Q.inf=Q.inf-diag(sum(Q.inf,2));
Q.QC=Q.QC-diag(sum(Q.QC,2));
Q.prog=Q.prog-diag(sum(Q.prog,2));

HHconfig.dataI = dataI;

return