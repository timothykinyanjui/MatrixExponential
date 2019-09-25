function kl = kl_Div(P,Q)
%
% Function calculates the KL divergence
% Written by Tim Kinyanjui
% University of Manchester 6th May 2015
% Input:
%   P: True distribution
%   Q: Theory/Appoximation
%
% K-L Divergence of Q from P measures the information gained when on
% revises ones beliefs from the prior (Q) to the posterior belief (P)

% Store the values
% kl_1 = zeros(length(P),1);
% 
% parfor i = 1:length(P)
%     
%     % If any is zero, store zero
%     if P(i)==0
%         
%         kl_1(i) = 0;
%         
%     elseif P(i)<0 || Q(i)<0 || P(i)>1 || Q(i)>1
%         
%         warning('Check that output is a probability vector')
%         kl_1(i) = 9556.21; % Check this
%         
%     elseif P(i)>0 && P(i)<=1 && Q(i)>=0 && Q(i)<=1
%         
%         kl_1(i) = P(i)*(log(P(i))-log(Q(i)));
%         
%     end
%     
% end
% 
% % Sum up
% kl = sum(kl_1);

% Do least squares
kl = sqrt(sum((P-Q).^2));

% End of function
return