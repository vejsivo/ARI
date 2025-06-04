%POLDEMOMIXSENS Demo for mixed sensitivity optimization with Polynomial Toolbox V2.5.
%
% Author(s):    Z. Hurak 
% Copyright (c) 2000 by PolyX, Ltd.
% $Revision: 1.0.1 $    $Date: 01-05-2001 $
%-----------------------------------------------------------------------------------------
clc
global PGLOBAL;
eval('PGLOBAL.FORMAT;', 'painit;');

format compact
echo on
%***************************************************************
%*                                                             *
%*   DEMO:  H_INFINITY MIXED SENSITIVITY OPTIMIZATION          *
%*                                                             *
%***************************************************************
 
pause %strike any key to continue

% A magnetic memory read/write head assembly can be well modelled 
% at low frequencies as
%
%           1
% G0(s) = -----
%          s^2
%
% The arm supporting the read/write facility head has some lightly 
% damped flexibility with uncertain resonant frequency. The more
% accurate model is
%
%                   s
%                B ---- + 1
%         30        w0
% G(s) = ----- ----------------------- 
%         s^2     s^2        s
%                ------ + B ---- + 1
%                 w0^2       w0
%
%
% We are expected to find a feedback compensator that guarantees closed
% loop stability for damping ratio B = 0.03 - 0.3 and resonant frequency 
% w0 >= 1 rad/sec.
%
% Moreover, it is required that the rise time of the closed loop system is 
% less than 1 sec and the overshoot is less than 5 %. The closed loop system
% should have zero steady state errors to ramps.
%
%---------------------------------------------------------------------------
pause %strike any key to continue...

% Both the uncertainty in the damping and the uncertainty in the resonant 
% frequency affect the high frequency plant characteristics. We model the as
% perturbation in numerator:
%          
%                       s^2
%                    - ------
%    N(s)-N0(s)         w0^2
%   ----------- = -------------------- = W(s)
%       N0(s)      s^2        s
%                 ------ + B ---- + 1
%                  w0^2       w0	

pause %strike any key to continue...

% We plot the magnitude of the frequency response of W for the worst-case 
% values of B and w0:

B = 0.03; w0 = 1;

numW = -s^2/w0^2;
denW = s^2/w0 + B*s/w0 + 1;

omega = logspace(-1,1);
magnW = bode(pol2mat(numW),pol2mat(denW),omega);
loglog(omega,magnW)
grid on, ylabel('|W(jw)|'), xlabel('w   [rad/sec]') 
title('Model of the numerator perturbation')

pause %strike any key to continue...

% Requirements on control ('translated'):
% ----------------------------------------
% 
% 1) Closed loop bandwidth 1 rad/sec
%
% 2) Small sensitivity at low frequencies, slope 2 decades/decade (40 db/dec) 
%    to guarantee type-2 control. 
%
% 3) Sufficiently fast decrease of the complementary sensitivity function 
%    above the frequency 1 rad/sec.
%
% 4) Relative damping of dominant poles 0.7, natural frequency 1 rad/sec.

pause %strike any key to continue

% The nominal plant is described:

num0 = 30; den0 = s^2;

% The two dominant poles are to be placed at s = -1/sqrt(2) -+ j*1/sqrt(2)

m = (s - 0.5*sqrt(2)*(-1+j))*(s - 0.5*sqrt(2)*(-1-j));

% The shaping function V(s) is then

numV = m;
denV = den0;

% To shape the sensitivity function S(s) at low frequencies we still have to
% determine the weighting function W1(s). Let's try W1(s)=1:

numW1 = 1;
denW1 = 1;

omega = logspace(-1,1);
magnInvVW1 = bode(pol2mat(denV*denW1),pol2mat(numV*numW1),omega);
loglog(omega,magnInvVW1)
grid on, ylabel('|1/(V(jw)*W1(jw))|'), xlabel('w   [rad/sec]') 
title('Approximation of the sensitivity function at low frequencies')

pause %strike any key to continue

% To shape the complementary sensitivity function T(s) at high frequencies, we must
% determine the weighting function W2(s). To consider it in the form W2(s) = c*(r*s+1).
% Experimenting a bit we choose c = 0.1. We chose r = 1 to guarantee high 
% frequency roll-off of T(s) of 3 decades/decade (60 dB/dec), the input sensitivity 
% function U(s) is then 20 dB/dec. Note, that this makes the weighting term nonproper:

c = 1;
r = 1;

numW2 = c*(r*s+1);
denW2 = 1;

% We compute the feedback compensator minimizing the infinity norm of the so-called
% general H_infinity system:
%
%        [ W1 S V ]
% H(s) = |        |
%        [ W2 U V ]

gmin = 1;                          %lower bound on the Hinf norm
gmax = 5;                          %upper bound on the Hinf norm
acc = 1e-4;                        %accuracy
tol = [1e-4,1e-8,5e-7,1e-8];       %tolerances (see the manual for more details)



[numC,denC] = mixeds(num0,m,den0,numW1,denW1,numW2,denW2,gmin,gmax,acc,tol)

pause %strike any key to continue

% the closed loop poles are

clpoles = roots(den0*denC + num0*numC)

% We see, that the dominant poles are at desired positions.

pause %strike any key to continue

omega = logspace(-2,2);
S = bode(pol2mat(den0*denC),pol2mat(den0*denC+num0*numC),omega);
T = bode(pol2mat(num0*numC),pol2mat(den0*denC+num0*numC),omega);

subplot(1,2,1), loglog(omega,S),
xlabel('\omega [rad/sec]');
ylabel('| S(j\omega) |');
axis([1e-2 1e2 1e-5 1e1]);
grid on;
subplot(1,2,2), loglog(omega,T)
xlabel('\omega [rad/sec]');
ylabel('| T(j\omega) |');
axis([1e-2 1e2 1e-5 1e1]);
grid on;
echo off

%end .. poldemomixsens
