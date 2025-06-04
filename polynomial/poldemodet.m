%POLDEMODET Comparision of numerical and symbolic computation of polynomial matrix determinant.
%
% Author(s):    M. Sebek, Z. Hurak 
% Copyright (c) 2000 by PolyX, Ltd.
% $Revision: 1.0 $    $Date: 12-13-2000 $
%-----------------------------------------------------------------------------------------
clc
global PGLOBAL;
eval('PGLOBAL.FORMAT;', 'painit;');

format compact

eval('syms tempVar; clear tempVar;',...
   'error(''Symbolic Math Toolbox required.'');');

echo on
%***************************************************************
%*                                                             *
%*   Comparision of numerical and symbolical computation       *
%*   of determinant of polynomial matrix.                      *
%*                                                             *
%***************************************************************
 
pause %strike any key to continue
 
% First, let's create some polynomial matrix. 
% We can use the PRAND function to generate a random square polynomial 
% matrix of dimension 6 x 6 and of degree 6.
  
pause %strike any key to continue
   
P=prand(6,6)
   
pause %strike any key to continue
 
% Now, lets create the same matrix in Symbolic Math Toolbox format.
% Towards this end, use the SYM function.

pause %strike any key to continue
    
Ps=sym(P)
    
pause %strike any key to continue
 
% We compute the determinant of the polynomial matrix using Symbolic 
% Math Toolbox and measure the time it takes:
 
tic,ps=det(Ps);toc
 
pause %strike any key to continue
 
% We also compute the determinant of the polynomial matrix numerically 
% using the Polynomial Toolbox:
 
tic,p=det(P);toc
   
pause %strike any key to continue
 
%The difference between the two results is then
 
dp = p-ps

echo off

%end .. poldemodet
