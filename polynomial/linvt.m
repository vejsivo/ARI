function Q = linvt(P, a, b);
%LINVT  Linear transformation of the variable of a polynomial matrix
%
% Given the polynomial matrix P with variable VAR and the real 
% numbers A and B the commmand
%    Q = LINVT(P,A,B)
% computes the polynomial matrix Q such that 
%	 Q(VAR) = P(A*VAR+B)
% If the third input is missing then B is assumed to be zero.
%
% This macro exists only fo completeness.
% See also POL/LINVT.

%       Author(s): M. Hromcik, M. Sebek 3-9-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 03-Sep-1998 10:28:34   $

if nargin == 0,
   error('Not enough input arguments.');
elseif nargin == 1,
  a = 1; b = 0;
elseif nargin == 2,
  b = 0;
end;
    
if ~isa(a, 'double') | ~isa(b, 'double') | ...
    any( size(a) > 1 ) | any( size (b) >1 ),  
  error('Invalid 2nd or 3rd argument; must be scalars.');
end;

if ~isa(P,'double') | ndims(P)>2,
   error('Invalid 1st argument.');
end;

Q = P;

%end .. linvt
