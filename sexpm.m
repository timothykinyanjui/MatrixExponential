function [expmA,s,k,m] = sexpm(A,shift,frac)
% [expmA,s,k,m] = sexpm(A,shift,frac)
% uses efficient subdiagonal Pade for expm(A), dense matices. 
% shift is estimate for rightmost eigval (optional)
%
% DISCLAIMER: this code is designed for matrices A with the properties 
% described in the reference below. In particular, it may produce rubbish results if the
% rightmost eigenvalues contain largely varying imaginary parts. 
%
%
% frac: character, optional input
% frac ='p': partial fraction form (default)
% frac ='v': product form (gets slightly better accuracy when ||A||<=O(1))
% frac ='r': use fixed parameters s=4 and (4,5) Pade for ||A||>1
%
% outputs: expm(A), scaling s, Pade degree (k,m)
%
%
%
% For details of the algorithm, see 
% S. Guettel and Y. Nakatsukasa, 
% Scaled and squared subdiagonal Pade
% approximation for the matrix exponential, 2015, preprint.

n = length(A);

if exist('shift','var')==0 
    shift = max(real(eig(A)));        
elseif isempty(shift);
    shift = max(real(eig(A)));            
end

if exist('frac')==0;
    frac='p'; % partial fraction
    %frac = 'v'; % product
end

A = A-shift*eye(n); % shift so that largest real part is (about) 0

% choose scaling/degree parameters to efficiency+stability
nrm = normest(A,2e-1);
if nrm>1
    if sum(frac=='r')>0, 
    s = 4; k = 4; m = k+1;   % fixed
    else
if nrm<200,      s = 4; k = 5; m = k-1;     
elseif  nrm<1e4, s = 4; k = 4; m = k+1;
elseif  nrm<1e6, s = 4; k = 3; m = k+1;
elseif  nrm<1e9, s = 3; k = 3; m = k+1;    
elseif  nrm<1e11,s = 2; k = 3; m = k+1;        
elseif  nrm<1e12,s = 2; k = 2; m = k+1;            
elseif  nrm<1e14,s = 2; k = 1; m = k+1;                
else s = 1; k = 1; m = k+1;                    
end
    end
else % nrm<1
if nrm>.5,      s = 4; k = 4; m = k-1;
elseif nrm>.3,  s = 3; k = 4; m = k-1;    
elseif nrm>.15, s = 2; k = 4; m = k-1;    
elseif nrm>.07, s = 1; k = 4; m = k-1;
elseif nrm>.01, s = 0; k = 4; m = k-1;    
elseif nrm>3e-4, s = 0; k = 3; m = k-1;
elseif nrm>1e-5,s = 0; k = 3; m = 0;
elseif nrm>1e-8,s = 0; k = 2; m = 0;
else s = 0; k = 1; m = 0;    % exp(A) = I+A to eps!
    expmA = exp(shift)*(eye(n)+A);
    return
end
end

% start method
As=A/(2^s); % scaled matrix

if frac=='p' | m==0 % partial fraction, default
[r,p,remainterm] = getcoeffs(k,m);

expmA = zeros(size(A)); % initialize
if isreal(A)
for ii=1:length(r)/2 % use partial fraction to get R(As)
RR = r(2*ii)*inv(As-p(2*ii)*eye(n));
expmA = expmA+RR+conj(RR);
end
if mod(length(p),2)==1
RR = r(end)*inv(As-p(end)*eye(n));    
expmA = expmA+RR;
end
else % A complex

for ii=1:length(r) % use partial fraction to get R(As)
RR = r(ii)*inv(As-p(ii)*eye(n));
expmA = expmA+RR;
end    
end

for ii=1:length(remainterm)
    RR = remainterm(end-ii+1)*As^(ii-1);
    expmA=expmA+RR;
end

else % product form        
    [rootp,rootq,mult] = getcoeffsproduct(k,m);
     % get coeffs for product form    
    expmA = eye(n);
    minlen = min(length(rootp),length(rootq));
    for ii=1:minlen
        expmA = (As-rootp(ii)*eye(n))*expmA;
        expmA = inv(As-rootq(ii)*eye(n))*expmA;
    end    
    % extra enumerator
    for ii = minlen+1:length(rootp)
        expmA = (As-rootp(ii)*eye(n))*expmA;
    end    
    % extra denominator
    for ii = minlen+1:length(rootq)
        expmA = inv(As-rootq(ii)*eye(n))*expmA;
    end
    expmA = expmA*mult;

end

for ii = 1:s % final squaring
expmA = expmA^2;
end
expmA = expmA*exp(shift);
end


function [r,q,remain] = getcoeffs(k,m) % table of coefficients for each case

if m == k+1; 
    remain = [];
    if k == 4
     r =  [ -1.582680186458572e+01 - 2.412564578224361e+01i;...
     -1.582680186458572e+01 + 2.412564578224361e+01i;...
      1.499984465975511e+02 + 6.804227952202417e+01i;...
      1.499984465975511e+02 - 6.804227952202417e+01i;     
     -2.733432894659307e+02                         ];...
     q = [   3.655694325463550e+00 + 6.543736899360086e+00i;...
      3.655694325463550e+00 - 6.543736899360086e+00i;...
      5.700953298671832e+00 + 3.210265600308496e+00i;...
      5.700953298671832e+00 - 3.210265600308496e+00i;...
        6.286704751729261e+00                        ];
    elseif k==3
  r = [-1.130153999597152e+01 + 1.247167585025031e+01i;...
     -1.130153999597152e+01 - 1.247167585025031e+01i;...
      1.330153999597152e+01 - 6.007173273704750e+01i;...
      1.330153999597152e+01 + 6.007173273704750e+01i];
  
   q=[3.212806896871536e+00 + 4.773087433276636e+00i;...
      3.212806896871536e+00 - 4.773087433276636e+00i;...
      4.787193103128464e+00 + 1.567476416895212e+00i;...
      4.787193103128464e+00 - 1.567476416895212e+00i];
  
    elseif k==2
   r=[7.648749087422928e+00 + 4.171640244747463e+00i;...
      7.648749087422928e+00 - 4.171640244747463e+00i;...
     -1.829749817484586e+01                         ];

      q = [2.681082873627756e+00 + 3.050430199247411e+00i;...
      2.681082873627756e+00 - 3.050430199247411e+00i;...
      3.637834252744491e+00           ];
    elseif k==1
    r = [ 1.000000000000000e+00 - 3.535533905932738e+00i;...
      1.000000000000000e+00 + 3.535533905932738e+00i];

    q = [2.000000000000000e+00 + 1.414213562373095e+00i;...
      2.000000000000000e+00 - 1.414213562373095e+00i];        
    end
    return
end
if m==k-1,
    if k==5
r = [     -1.423367961376821e+02 - 1.385465094833037e+01i;...
     -1.423367961376821e+02 + 1.385465094833037e+01i;...
      2.647367961376822e+02 - 4.814394493714596e+02i;...
      2.647367961376822e+02 + 4.814394493714596e+02i];
q = [      5.203941240131764e+00 + 5.805856841805367e+00i;...
      5.203941240131764e+00 - 5.805856841805367e+00i;...
      6.796058759868242e+00 + 1.886649260140217e+00i;...
      6.796058759868242e+00 - 1.886649260140217e+00i];
  
remain =   [2.000000000000000e-01     9.8000000000000e+00];
    elseif k==4        
r = [    2.484269593165883e+01 + 7.460342395992306e+01i;...
      2.484269593165883e+01 - 7.460342395992306e+01i;...
     -1.734353918633177e+02                         ];
q = [      4.675757014491557e+00 + 3.913489560603711e+00i;...
      4.675757014491557e+00 - 3.913489560603711e+00i;...
      5.648485971016893e+00                         ];
remain =[    -2.500000000000000e-01    -7.750000000000000e+00        ];
    elseif k==3
r=  [    2.533333333333333e+01 - 2.733333333333333e+01i;...
      2.533333333333333e+01 + 2.733333333333333e+01i];
q = [      4.00000000000000e+00 + 2.00000000000000e+00i;...
      4.00000000000000e+00 - 2.00000000000000e+00i];
remain =[     3.333333333333333e-01     5.666666666666667e+00        ];
    elseif k==2        
r =  -13.5     ;
q =    3       ;
remain =[ -0.5  -3.5];
    end    
end
if m==0
    r=[];q=[];    
    if k==3
 remain = [1/6            1/2            1              1           ];
    elseif k==2
 remain = [   1/2            1              1       ];        
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%

function [p,q,mult] = getcoeffsproduct(k,m)
if m == k+1; 
      mult = (-1)^(m)*m;
    if k == 4
        
p = [ -5.203941240131764e+00 + 5.805856841805367e+00i;...
     -5.203941240131764e+00 - 5.805856841805367e+00i;...
     -6.796058759868242e+00 + 1.886649260140217e+00i;...
     -6.796058759868242e+00 - 1.886649260140217e+00i];

q = [ 3.655694325463550e+00 + 6.543736899360086e+00i;...
      3.655694325463550e+00 - 6.543736899360086e+00i;...
      6.286704751729261e+00                         ;...
      5.700953298671832e+00 + 3.210265600308496e+00i;...
      5.700953298671832e+00 - 3.210265600308496e+00i];
  
    elseif k==3
p=[     -4.675757014491557e+00 + 3.913489560603711e+00i;...
     -4.675757014491557e+00 - 3.913489560603711e+00i;...
     -5.648485971016893e+00                         ];

q=[      3.212806896871536e+00 + 4.773087433276636e+00i;...
      3.212806896871536e+00 - 4.773087433276636e+00i;...
      4.787193103128464e+00 + 1.567476416895212e+00i;...
      4.787193103128464e+00 - 1.567476416895212e+00i    ];    
    elseif k==2
    p = [-4.00000000000000e+00 + 2.00000000000000e+00i;...
     -4.00000000000000e+00 - 2.00000000000000e+00i];
     q = [2.681082873627756e+00 + 3.050430199247411e+00i;...
      2.681082873627756e+00 - 3.050430199247411e+00i;...
      3.637834252744491e+00                         ];        
    elseif k==1
   p = -3;
   q = [2.000000000000000e+00 + 1.414213562373095e+00i;...
      2.000000000000000e+00 - 1.414213562373095e+00i      ];  
    end
    return
end
if m==k-1,
  mult = (-1)^(m)/k    ;
    if k==5
p=  [   -3.655694325463550e+00 + 6.543736899360086e+00i;...
     -3.655694325463550e+00 - 6.543736899360086e+00i;...
     -6.286704751729261e+00                         ;...
     -5.700953298671832e+00 + 3.210265600308496e+00i;...
     -5.700953298671832e+00 - 3.210265600308496e+00i];

q=  [    5.203941240131764e+00 + 5.805856841805367e+00i;...
      5.203941240131764e+00 - 5.805856841805367e+00i;...
      6.796058759868242e+00 + 1.886649260140217e+00i;...
      6.796058759868242e+00 - 1.886649260140217e+00i    ];    
    elseif k==4        
p= [   -3.212806896871536e+00 + 4.773087433276636e+00i;...
     -3.212806896871536e+00 - 4.773087433276636e+00i;...
     -4.787193103128464e+00 + 1.567476416895212e+00i;...
     -4.787193103128464e+00 - 1.567476416895212e+00i];
q= [     4.675757014491557e+00 + 3.913489560603711e+00i;...
      4.675757014491557e+00 - 3.913489560603711e+00i;...
      5.648485971016893e+00                             ];      
    elseif k==3
    p= [ -2.681082873627756e+00 + 3.050430199247411e+00i;...
     -2.681082873627756e+00 - 3.050430199247411e+00i;...
     -3.637834252744491e+00                         ];
    q =[ 4.00000000000000e+00 + 2.00000000000000e+00i;...
      4.00000000000000e+00 - 2.000000000000001e+00i];        
    elseif k==2        
   p = [-2.000000000000000e+00 + 1.414213562373095e+00i;...
     -2.000000000000000e+00 - 1.414213562373095e+00i];
   q= 3        ;
    end    
end
end
