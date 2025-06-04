function p = minsens(N,D)
%MINSENS  Minimum peak sensitivity
%
% The function
%    p = minsens(N,D)
% computes the minimum peak value p of the sensitivity and
% complementary sensitivity functions for the SISO plant 
% with transfer function P = N/D that may be achieved by
% feedback.
%
% WARNING: This macro is part of the demo "Achievable feedback 
% performance." See the Manual for an explanation. No error checks 
% have been built in.

% Reference:
% Kwakernaak, H. (1985). "Minimax frequency domain performance 
% and robustness optimization of linear feedback systems." 
% IEEE Transactions on Automatic Control, vol. AC-30, pp. 994-1004.

% H. Kwakernaak, January, 1999
% Copyright PolyX Ltd, 1999

if nargin~=2
   disp('minsens: This macro is part of the demo "Achievable feedback') 
   disp('         performance." See the Manual for an explanation.') 
   return
end

% Compute the polynomials Nplus and Dplus whose roots are 
% the roots of N and D, respectively, with positive real parts

rootsN = roots(N); rootsNplus = rootsN(find(real(rootsN)>0));
Nplus = mat2pol(poly(rootsNplus));
rootsD = roots(D); rootsDplus = rootsD(find(real(rootsD)>0));
Dplus = mat2pol(poly(rootsDplus));
n = deg(Nplus); d = deg(Dplus); 


% Check wether p = 0 or p = 1

if n == 0
   p = 0; return
elseif d == 0
   p = 1; return
end


% Solve the generalized eigenvalue problem

A = [ sylv(Dplus,'col',n-1) sylv(Nplus,'col',d-1) ];
J = 1;
for i = 2:n
   J(i,i) = -1*J(i-1,i-1);
end
B = [ sylv(Dplus','col',n-1)*J zeros(n+d,d) ];
p = 1/min(abs(eig(A,B)));

%end .. minsens
