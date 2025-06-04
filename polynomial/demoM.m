% demoM  Script file for the demo "Polynomial solution 
%        of the SISO mixed sensitivity H-infinity problem." 
%        See the Manual for an explanation. Assign a value 
%        to the scalar gamma before running the script.

% Huibert Kwakernaak, January, 1999
% Modified by S. Pejchova, June 9, 1999.
% Copyright PolyX Ltd, 1999

disp('demoM  Script file for the demo "Polynomial solution') 
disp('       of the SISO mixed sensitivity H-infinity problem."') 
disp('       See the Manual for an explanation.')
disp('       Hit any key to continue')
pause

% Check the existence of gamma
if exist('gamma')~=1,
   error('Assign a value to gamma before running the script.');
end

% Define the data
n = 1; d = s^2; m = s^2+s*sqrt(2)+1; 
c = 1; r = 1;  a1 = 1; b1 = 1; a2 = c*(1+r*s); b2 = 1;  

% Define the various polynomial matrices
D1 = [b1 0; 0 b2; 0 0];
D2 = [a1  ; 0   ; d  ];
N1 = [ 0;  0; -m ];
N2 = [ 0; a2; -n ];  

% Left-to-right conversion
Q = diag([b1 b2 m]);
[Del,Lam] = lmf2rmf([N2 D2],Q);  

% Spectral factorization
Jgamma = eye(3); Jgamma(3,3)= -gamma^2;
DelDel = Del'*Jgamma*Del;  
% [Gam,J] = spf(DelDel); 
[Gam,J] = spf(DelDel,'nnc');

% Computation of the compensator
xy = rmf2lmf([1 0]*Gam,Lam);
x = xy(1,1); y = xy(1,2);

% Computation of the closed-loop characteristic polynomial 
% and closed-loop poles
phi = d*x+n*y;
clpoles = roots(phi)  

