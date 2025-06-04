% demoB  Script file for the demo "Control of a batch process"
%        See the Manual for an explanation.

% Huibert Kwakernaak, January, 1999
% Modified by S. Pejchova, June 9, 1999.
% Copyright PolyX Ltd, 1999

disp(' ');
disp('demoB: Script file for the demo "Control of a batch process."')
disp('       See the Manual for an explanation --- Hit any key to continue.');
pause

% Save workspace variables
save demoB;
ws=warning; warning off;

% Define the data
L = [-77*z^-3 0.33*z^-3; 0 0.03*z^-3; 0 0.1*z^-3];
R = eye(2,2);
theta = diag([1-0.6*z^-1 1 1-0.55*z^-1]);
Phi = diag([1 1-0.5*z^-1 1]);
d = 1;
Q1 = diag([0.01 8 2.25]);
Nabla = (1-z^-1);

% Spectral factorization
Gamma = spf(L'*Q1*L);

% Solution of the two-sided equation
n = 3;
A = z^-n*Gamma';
B = Nabla^d*Phi;
C = L'*Q1*theta;
[T,U] = axybc(A,B,C);

% Compute the controller H = phi/theta
phi = R/Gamma*T;

% Compute the sensitivity matrix S = psi/theta
psi = theta-L/Gamma*T;

% Compute and plot the disturbance impulse response matrix r
% Apply long division to psi/thet
n = 10;
[q,r] = longrdiv(psi,theta,n);
% Plot
figure; clf
k = length(r);
for i = 1:k
   for j = 1:k
      subplot(k,k,(i-1)*k+j)
      plot(0:n,r{:}(i,j))
      grid on; axis([0 n -1.5 1.5])
   end
end

% Compute and plot the disturbance step response matrix s
% Apply long division to (psi/theta)*1/(1-z^-1)
n = 10;
[q,s] = longrdiv(psi,theta*(1-z^-1),n);
warning(ws);

% Plot
figure; clf
k = length(s);
for i = 1:k
   for j = 1:k
      subplot(k,k,(i-1)*k+j)
      plot(0:n,s{:}(i,j))
      grid on; axis([0 n -1.5 1.5])
   end
end

%Refresh the workspace variables
clear A B C Gamma L Nabla Phi Q1 R T U d i j k n phi psi q r s theta ws
load demoB;
delete demoB.mat;
