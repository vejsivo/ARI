function [D,T] = tdeg(A,string)
%TDEG    Trailing degree of constant
%
% D = TDEG(A)       default, the same as TDEG(A,'mat').
% D = TDEG(A,'mat') returns the trailing degree of matrix A.
% D = TDEG(A,'ent') returns the matrix of trailing degrees of the A entries.
% D = TDEG(A,'row') returns the column vector of row trailing degrees of A.
% D = TDEG(A,'col') returns the row vector of column trailing degrees of A.
% D = TDEG(A,'dia') for symmetric matrix A,
%                    returns the vector of half diagonal trailing degrees of A.
%
% The second output argument T in all cases returns corresponding
% matrix or vector of trailing coefficients, the same as first output argument
% in function TCOEF.
%
% See also TCOEF, POL/TDEG, POL/TCOEF, TSP/TDEG, TSP/TCOEF.

%          Author: J. Jezek  13-10-99
%          Copyright (c) 1999 by Polyx, Ltd.
%          $ Revision 3.0 $   $ 15-Jun-2000 $

% Effect on other properties:
% D and T are standard Matlab matrices.

if nargin<1,
   error('Not enough input arguments.');
elseif nargin==1,
   string = 'mat';
elseif ~ischar(string),
   error('Invalid 2nd argument.');
end;

eval('[D,T] = tdeg(pol(A),string);','error(peel(lasterr));');

%end .. tdeg
