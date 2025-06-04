%POLDEMOROBPAR Stability analysis of systems with parametric uncertainties.
%
% Author(s):    M. Sebek, Z. Hurak 
% Copyright (c) 2000 by PolyX, Ltd.
% $Revision: 1.0 $    $Date: 12-13-2000 $
%-----------------------------------------------------------------------------------------

global PGLOBAL;
eval('PGLOBAL.FORMAT;', 'painit;');
clc
echo on
%**************************************************************************
%*                                                                        *
%* EXAMPLE 1: ROBUST STABILITY INTERVAL FOR ONE UNCERTAIN PARAMETER       *
%*                                                                        *
%**************************************************************************
 
% We have a characteristic polynomial whose coefficients depend affinely
% on one uncertain parameter. Then such a polynomial can be written:
%
% p(s,q)=p0(s)+q*p1(s) 
%
% where p0(s) is a stable polynomial,  q is the uncertain parameter. 
% Now, find the maximum stability interval [qmin, qmax].

pause %strike any key to continue

% First, enter the two polynomials p0(s) and p1(s)

p0 = 3 + 10*s + 12*s^2 + 6*s^3 + s^4; 
p1 = s + s^3;  

pause %strike any key to continue

% stability of p0(s) is easily tested:

isstable(p0)

pause %strike any key to continue

% Finding the stability interval is accomplished with the command STABINT:

[qmin,qmax]=stabint(p0,p1)  
 
pause %strike any key to continue

% We can check this result by plotting the root locus of a (artificial) plant p1/p0
% and feedback gain q. 
% This requires that the Control Toolbox 4.2 or higher be installed on your computer:
   
pause %strike any key to continue
   
rlocus(ss(p1,p0),qmin:.1:100)
   
pause %strike any key to continue
   
close
clc
%****************************************************************************
%*                                                                          *
%* EXAMPLE 2: CONTINUOUS-TIME INTERVAL POLYNOMIALS BY KHARITONOV THEOREM    *
%*                                                                          *
%****************************************************************************

% Consider the continous-time characteristic polynomial, whose coefficients
% are only know to lie within some intervals:
%
% p(s,q)=[0.45,0.55] + [1.95,2.05]*s + [2.95,3.05]*s^3 + [3.95,4.05]*s^4 + [3.95,4.05]*s^5 + s^6
%
% Check if it is robustly stable.

pause %strike any key to continue

% Enter the interval polynomial via two 'lumped' polynomials:

pminus = 0.45+1.95*s+2.95*s^2+5.95*s^3+3.95*s^4+3.95*s^5+s^6;
pplus = 0.55+2.05*s+3.05*s^2+6.05*s^3+4.05*s^4+4.05*s^5+s^6; 

pause %strike any key to continue

% Theorem of Kharitonov says that it is sufficient to test for the stability 
% of the four distinguished polynomials. To obtain these so-called Kharitonov 
% polynomials we use the KHARIT function:

[stability,K1,K2,K3,K4]=kharit(pminus,pplus)

pause %strike any key to continue
close
clc
%****************************************************************************
%*                                                                          *
%* EXAMPLE 3: CONTINUOUS-TIME INTERVAL POLYNOMIALS BY GRAPHICAL METHODS     *
%*                                                                          *
%****************************************************************************

% For the same data as in the previous example, use the graphical 
% method of Zero Exclusion Principle to test the robust stability 
% of the uncertain polynomial.

pause %strike any key to continue

% First, check if there is at least one stable member in the polynomial family:

isstable(pminus)

pause %strike any key to continue

% Plot the value sets (here Kharitonov rectangles) on imaginary axis.
% Let's make the initial guess at he frequency range OMEGA between 0 and 1:

khplot(pminus,pplus,0:.002:1)

pause %strike any key to continue
close
clc
%**************************************************************************
%*                                                                        *
%* EXAMPLE 4: DISCRETE-TIME INTERVAL POLYNOMIALS BY GRAPHICAL METHODS     *
%*                                                                        *
%**************************************************************************

% Even though there is no discrete-time counterpart to Kharitonov theorem,
% we can solve these problems efficiently in a general polytopic framework

% The polytope of polynomials is defined as 
% 
% p(z,q)=p0(z)+ q1*p1(z) + ... + qN*pN(z)
%
% where p0(z) must be stable, the bounds on parameters q1,...,qN are given 
% by the matrix 
%
% Qbounds = [q1min  q2max
%            q2min  q2max
%              .      .
%              .      .
%            qNmin  qNmax]
%
% We want to test, if p(z,q) is stable for all admissible values of q1,...,qN.  

pause %strike any key to continue


% Let's consider a discrete-time interval polynomial given by:
%
% p(z,q) = [10,20] + [20,30]*z + [128,138]*z^2 + [260,270]*z^3 + [168,178]*z^4

pause %strike any key to continue

% It can be described by the extreme polynomials:

p0 = 10 + 20*z + 128*z^2 + 260*z^3 + 168*z^4;
p1=1;
p2=z;
p3=z^2;
p4=z^3;

% and the bounds on the uncertain parameters are given by

Qbounds=[0 10;0 10;0 10;0 10];  

pause %strike any key to continue

% Now, we can easily plot the value sets (octagons) for generalized frequencies 
% on unit circle using the function PTOPPLOT:

ptopplot(p0,p1,p2,p3,p4,Qbounds,exp(j*(0.3:0.002:0.7)*2*pi))

% we see that zero is excluded from the value set, so we conclude that the uncertain 
% polynomial is robustly stable.

pause %strike any key
close
clc
%**************************************************************************
%*                                                                        *
%* EXAMPLE 5: GENERAL POLYTOPE OF POLYNOMIALS                             *
%*                                                                        *
%**************************************************************************

% As in the previous example, we can define even more general polynomial families 
% by entering the set of its extreme members, so-called generators of the polytope.
% This is especially useful, when the uncertain parameters enter the coefficients 
% of the uncertain polynomial affinely. 

pause %strike any key to continue

% Consider the following uncertain polynomial 
%
% p(s,q) = 3+q1-q2+2*q3 + (3+3*q1+q2+q3)*s + (3+3*q1-3*q2+q3)*s^2 + (1+2*q1-q2+2*q3)*s^3
% 
% with       -0.245 <= q1 <= 0.245
%            -0.245 <= q2 <= 0.245
%            -0.245 <= q3 <= 0.245
%
% Is the polynomial p(s,q) robustly stable?

pause %strike any key to continue

% First, enter the data:

p0=pol([3 3 3 1],3);
p1=pol([1 3 3 2],3);
p2=pol([-1 1 -3 -1],3);
p3=pol([2 1 1 2],3);
Qbounds=[-0.245 0.245;-0.245 0.245;-0.245 0.245]; 

pause %strike any key to continue

% Plot the value sets (octagons) on imaginary axis 
% with the initial guess at the range of frequencies from 0 to 1.5

ptopplot(p0,p1,p2,p3,Qbounds,j*(0:1.5/50:1.5))

% we see that zero is excluded from the value set (it can be verified that it holds
% even if we increase the frequency range, so we conclude that the uncertain 
% polynomial is robustly stable.

pause %strike any key to continue
close
clc
%*****************************************************************************
%*                                                                           *
%* EXAMPLE 6: UNCERTAIN POLYNOMIALS WITH MULTILINEAR UNCERTAINTY STRUCTURE   *
%*                                                                           *
%*****************************************************************************

% If the uncertain parameters enter the coefficients in multilinear way, 
% the nice theory from the previous examples doesn't apply. There are essentially
% two approaches to robust stability analysis of multilinear systems.
%
% 1) Overbounding by the polytope of polynomials and using the tools from the previous
%    examples. This however, introduces some conservatism into analysis. 
%
% 2) Gridding the set of uncertain parameters and plotting the value set polygons for 
%    a range of frequencies. Visually checking if the zero is excluded from the value 
%    set, we can make nonconservative conclusion about robust stability. Since this is
%    a 'brute force' approach it will be very time comsuming for high order systems. 

pause %strike any key to continue

% Consider the uncertain polynomial with multilinear uncertainty structure described by
%
% p(s,q) = p0(s)+(q1+q2)*p2(s)+q1*q2*p1(s) 
%
% where

p0=1+s+s^2;
p1=s;
p2=1;

% and the uncertain parameters range from 0 to 1:

q1=0:.01:1;
q2=0:.01:1;

pause %strike any key to continue

% Computation of the value sets for a range of frequencies (chosen from 0 to 2)
% is accomplished by the VSET and VSETPLOT functions:

V=vset(q1,q2,'p0+(q1+q2)*p2+q1*q2*p1',p0,p1,p2,j*(0:.1:2));vsetplot(V)

pause %strike any key to continue
close
clc
%*****************************************************************************
%*                                                                           *
%* EXAMPLE 7: UNCERTAIN POLYNOMIALS WITH POLYNOMIAL UNCERTAINTY STRUCTURE    *
%*                                                                           *
%*****************************************************************************

% The polynomial family exhibits a polynomial uncertainty structure 
% when the uncertain parameters enter in higher powers into the coefficients 
% of the polynomial. We can use the same tool as for multilinear case - the two 
% functions VSET and VSETPLOT accomplish the computation and plotting of the value
% sets by simple gridding of the set of uncertain parameters. 

pause %strike any key to continue

% Consider the uncertain polynomial
%
% p(s,q) = p0(s)+(q1^3+q1)*p1(s)+(q2^3-q2^2)*p2(s)+q1*q2*p3(s)
%
% with the two uncertain parameters q1 and q2 ranging from -1 to 1. Is the polynomial
% family robustly stable?

pause %strike any key to continue


uncertStruct='p0+(q1^3+q1)*p1+(q2^3-q2^2)*p2+q1*q2*p3'; 
p0=1+s+s^2;
p1=1+s;
p2=1-s;
p3=1-s+s^2;
q1=-1:.02:1; q2=-1:.02:1;

pause %strike any key to continue

% Computing and plotting the value sets on imaginary axis for frequencies from 0 to 4
% is easy with the VSET and VSETPLOT functions

V=vset(q1,q2,uncertStruct,p0,p1,p2,p3,j*(0:.4:4));vsetplot(V)

pause %strike any key to continue
close
clc
%*****************************************************************************
%*                                                                           *
%* EXAMPLE 8: UNCERTAIN POLYNOMIALS WITH GENERAL UNCERTAINTY STRUCTURE       *
%*                                                                           *
%*****************************************************************************

% The two functions for computing and plotting the values sets for a range of frequencies
% can be used for robust stability analysis of polynomial families with general uncertainty
% structures. Actually, the 'gridding' approach of the two functions will be the only
% solutino in such cases. Even though it might be rather time consuming, it is the 
% only tool available for analysis of such uncertain systems.

pause %strike any key to continue

uncertStruct='p0+sqrt(abs(q1))*p1+exp(q2)*p2+sin(q1*q2)*p3'; 
p0=1+s+s^2;
p1=1+s;
p2=1-s;
p3=1-s+s^2;
q1=-1:.02:1; q2=-1:.02:1;

pause %strike any key to continue

% Computing and plotting the value sets on imaginary axis for frequencies from 0 to 4
% is done by the VSET and VSETPLOT functions

V=vset(q1,q2,uncertStruct,p0,p1,p2,p3,j*(0:.4:4));vsetplot(V)

pause %strike any key to continue

% the uncertainty structure described in the string UNCERTSTRUCT is a very flexible way
% of entering whatever kind of dependence of coefficients on the uncertain parameters.
% 
% Let's ilustrate this by presenting a few more examples of uncertainty structures,
% using the same data.

pause %strike any key to continue

uncertStruct='p0+sqrt(abs(q1))*p1+sin(q2)*p2+cos(q1*q2)*p3'; 
V=vset(q1,q2,uncertStruct,p0,p1,p2,p3,j*(0:.4:4));vsetplot(V)

pause %strike any key to continue

uncertStruct='p0+q2*p1+sinh(q2)*p2+cosh(q1*q2)*p3'; 
V=vset(q1,q2,uncertStruct,p0,p1,p2,p3,j*(0:.4:4));vsetplot(V)

pause %strike any key to continue

uncertStruct='p0+airy(q1)*p1+airy(q2)*p2+(q1*q2)*p3'; 
V=vset(q1,q2,uncertStruct,p0,p1,p2,p3,j*(0:.4:4));vsetplot(V)

pause %strike any key to continue
close
clc
%*****************************************************************************
%*                                                                           *
%* EXAMPLE 9: SPHERICAL POLYNOMIAL FAMILIES                                  *
%*                                                                           *
%*****************************************************************************

% If the admissible set of uncertain parameters is defined by bounds on L2 norm,
% such polynomials are called spherical uncertain polynomials. SPHERPLOT function 
% is a tool for computing and plotting the value sets.

pause %strike any key to continue

% Consider the polynomial family given by
%
% p(s,q) = p0(s) + q1*s + q2*s + q3*s
%
% where p0(s) = 0.5 + s + 2*s^2 + 4*s^3and the weighted L2 norm of the vector 
% of uncertain parameters is less than 1. The weighting matrix for the norm is given by 
%      
% W = [2 0 0 0
%      0 5 0 0
%      0 0 3 0
%      0 0 0 1]
%
% Is the polynomial family robustly stable?

pause %strike any key to continue

% First enter the nominal polynomial, vector of diagonal entries of the weighting matrix,
% the bound on the L2 norm and the frequency range:

p0 = 0.5 + s + 2*s^2 + 4*s^3;
weight = [2 5 3 1]; 
bound = 1;
omega = 0:.025:1.5;

pause %strike any key to continue

% The value sets are obtained with the SPHERPLOT function:
spherplot(p0,omega,bound,weight);

pause %strike any key to continue
close
echo off

%end .. poldemorobpar
