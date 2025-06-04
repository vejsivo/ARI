function r = covf(A,C,n)
%COVF  Covariance function of an ARMA process
%
% The function call
%    R = COVF(A,C,N)
% computes and plots the covariance function r of the 
% discrete-time ARMA process y defined by
%    A(z)y(t) = C(z)e(t)
% with e standard white noise, up to the time shift N.
%
% WARNING: This macro is part of the demo "Computing the
% covariance function of an ARMA process." It has no input
% checks. Please see the Manual for an explanation.

% H. Kwakernaak, January, 1999
% Modified by S. Pejchova, June 9, 1999.
% Copyright by PolyX Ltd, 1999

if nargin~=3
   disp('covf: This macro is part of the demo "Computing the')
   disp('      covariance function of an ARMA process." See')    
   disp('      the Manual for an explanation.')
   return
end

% Solve the two-sided polynomial matrix equation
X = xaaxb(A,C*C');

% Apply long division to X\A
[Q,R] = longldiv(X,A,n);

% Construct the covariance function
r = R;
r.c(:,:,1)=r.c(:,:,1)+(r.c(:,:,1))';

% Plot the covariance function
figure; clf
k = length(C);
for i = 1:k
   for j = 1:k
      subplot(k,k,(i-1)*k+j)
      rij=r.c(i,j,:);
      plot(0:n,rij(:)','o')     
      grid on; xlabel('tau')
   end
end

%end .. @pol/covf
