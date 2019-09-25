function P = daMethod(h,II,M,T,order,PP0)

% Function solves using the Backward Euler method
% Written by Tim on 7th Aug 2015
% University of Manchester

% Set time vector
time = 0:h:T;

% Solve using Euler backward of order order
P(:,1) = PP0';
Msum = II - h*M;

if order >= 2
    for k = 2:order
        Msum = Msum + (((-1)^k)*(h^k)/prod(1:k))*M^k;
    end
end

for i = 2:length(time)
    P(:,i) = Msum\P(:,i-1);
end
P = P';

% End of function
return

% % Allows the modification of the step size to fit tolerance
% flag = true;
% 
% while flag
%     
%     % Set time vector
%     time = 0:h:T;
%     
%     % Solve using Euler backward of order order
%     P(:,1) = PP0';
%     Msum = II - h*M;
%     
%     if order >= 2
%         for k = 2:order
%             Msum = Msum + (((-1)^k)*(h^k)/prod(1:k))*M^k;
%         end
%     end
%     
%     for i = 2:length(time)
%         P(:,i) = Msum\P(:,i-1);
%     end
%     P = P';
%     
%     % Calculate to check sol (last 10 values) are within tolerance
%     nEnd = 5; % Last 10 points
%     lP = length(P(:,1));
%     parfor i = 0:(nEnd-1)
%         err(i+1,1) = max(abs(P((lP-i),:) - P((lP-(i+1)),:)));
%     end
%     
%     % Check tolerance
%     errTol = sum(err > tol);
%     
%     if errTol >= 1
%         fprintf('Reducing time step to %f to get tolerance of within %d\n',h,tol)
%         h = h/2;
%         flag = true;
%         clear P
%     else
%         flag = false;
%         disp('Done')
%     end
% 
% end