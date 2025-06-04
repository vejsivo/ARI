%POLDEMODEBE deadbeat compensator design
% 
% Author(s):    M. Sebek, Z. Hurak 
% Copyright (c) 2000 by PolyX, Ltd.
% $Revision: 1.0 $    $Date: 12-13-2000 $
%-----------------------------------------------------------------------------------------

global PGLOBAL;
eval('PGLOBAL.FORMAT;', 'painit;');
clc
format compact
echo on
%**************************************************************************
%*                                                                        *
%* EXAMPLE 1: DESIGN OF A DEADBEAT COMPENSATOR FOR A SISO PLANT           *
%*                                                                        *
%**************************************************************************

% Consider the SISO plant described by the transfer function 
%
% G(z) = num(z)/den(z)
% 
% where

num=1+z+z^2;
den=4+3*z+2*z^2+z^3;

pause %strike any key to continue

% The design of a deadbeat compensator is performed by the DEBE function:

[numc,denc]=debe(num,den)  

pause %strike any key to continue

% The closed-loop characteristic polynomial is then

cl=den*denc+num*numc

pause %strike any key to continue

% The step response of such system is obtained with the help of the Control Toolbox:

step(feedback(ss(num,den),ss(numc,denc)))

pause %strike any key to continue
close

% It is also possible to try to find other solutions of order 5
[nc,dc,e,f,degt]=debe(num,den)

pause %strike any key to continue

% The value of the last output parametr reveals that it is impossible to find 
% other proper deadbeat compensator of order 5. If we insist on finding such 
% a deadbeat compensators, we can solve the associated Diophantine equation directly:

[nc,dc,f,e] = axbyc(num,den,z^6)

% which describes the set of sixth order compensators parametrized by a constant T.

pause %strike any key to continue

% Alternatively, we can perform the same design in the backward-shift operator z^-1.

pause %strike any key to continue

[nneg,dneg]=reverse(num,den); symbol(nneg,'z^-1'); symbol(dneg,'z^-1'); nneg,dneg

[ncneg,dcneg]=debe(nneg,dneg)  

pause %strike any key

% The closed-loop characteristic polynomial is then

clneg=dneg*dcneg+nneg*ncneg

pause %strike any key
clc
%**************************************************************************
%*                                                                        *
%* EXAMPLE 2: DESIGN OF A DEADBEAT COMPENSATOR FOR A MIMO PLANT           *
%*                                                                        *
%**************************************************************************

% Consider the MIMO plant described by the right matrix fraction 
%
% G(z) = NUM(z)*DEN^-1(z)
% 
% where

NUM=[1-z z; 2-z 1],DEN=[1+2*z-z^2 -1+z+z^2; 2-z 2+3*z+2*z^2]

pause %strike any key to continue

% Design of a deadbeat compensator is performed by the DEBE function:

[NUMC,DENC]=debe(NUM,DEN,'r')   

pause %strike any key to continue

% The closed-loop characteristic polynomial is then

cl=det(DENC*DEN+NUMC*NUM)

% The distribution of poles over invariant polynomials can be seen from the Smith form:

smith(DENC*DEN+NUMC*NUM)

pause %strike any key to continue
echo off
format

%end .. poldemodebe
