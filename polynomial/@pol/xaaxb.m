function X = xaaxb(A,B,varargin)
%XAAXB  Symmetric polynomial equation solver
%
% The commmand
%    X = XAAXB(A,B) 
% solves the bilateral symmetric matrix polynomial equation
%    XA' + AX' = B
% where B is a para-Hermitian polynomial matrix, that is,
%    B = B' .
%
% For polynomials in variable 's' or 'p', see POL/CTRANSP.
% For polynomials in variable 'z' or 'z^-1', the symmetric B
% is a two-sided polynomial, see TSP/TSP, POL/CTRANSP, TSP/TRANSP.
% The argument B should be of class TSP, see TSP/AXXAB.
%
% For compatibility with the older version of the Polynomial Toolbox,
% for discrete-time polynomial, the symmetic argument B can be also
% of class POL (or convertible to that). In such a case, 
% the symmetry condition runs
%    B = BH' + BH 
% The degree offset in matrix B is evaluated upon cancelation of
% the leading and trailing zero matrix coefficients. If there are no
% zero coefficients then the degree offset is degB/2.
%
% The commmand
%    XAAXB(A,B,'syl') 
% solves the equation with the Sylvester matrix method. This is the 
% default method. It may be used with several modifiers:
%
% The commmand
%    XAAXB(A,B,'tri') 
% returns a solution with a lower-triangular row-wise leading coefficient 
% matrix (continuous-time case) or a lower-triangular absolute coefficient  
% matrix (discrete-time case).
%
% In the continuous-time case 
%    XAAXB(A,B,ROWDEG) 
% computes a solution X of row degrees ROWDEG. By default, 
%    ROWDEG(i) = MAX(DEG(B)-DEG(A),0)
% for all row indices i.
%
% In the discrete-time case the macro computes a solution X of degree
% DEG(X) = MAX(DEG(A), DEG(B)).
%    
% The commmand
%    XAAXB(A,B,'red') 
% solves the equation with polynomial reductions, a version of the Euclidean 
% division algorithm for polynomials. If A is stable and A^-1 B A'^-1 is biproper
% then the macro computes a solution with lower-triangular row-wise leading 
% coefficient matrix (continuous-time) and such that X^-1A is proper.
%
% If there is no solution of degree as specified above then all the entries in X 
% are set equal to NaN.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: TSP/XAAXB, POL/AXXAB, TSP/AXXAB.

%    Author: D. Henrion, September 10, 1998.
%    Copyright 1998 by Polyx, Ltd.

% Call AXXAB with tranposed arguments.

if nargin < 2,
 error('Not enough input arguments.');
elseif nargin < 3,
 varargin = [];
end;

eval('X = axxab(A'',B,varargin)'';','error(peel(lasterr));');

%end .. @pol/xaaxb


