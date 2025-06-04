function [T,D] = tcoef(A,string)
%TCOEF   Trailing coeficient of polynomial
%
% T = TCOEF(A)         default, the same as TCOEF(A,'mat')
% T = TCOEF(A,'mat')   returns the trailing coefficient matrix of A
% T = TCOEF(A,'ent')   returns the scalar trailing coefficients
%                        of the entries of A
% T = TCOEF(A,'row')   returns the row trailing coefficient matrix of A
% T = TCOEF(A,'col')   returns the column trailing coefficient matrix of A
% T = TCOEF(A,'dia')   for para-Hermitian polynomial matrix A,
%                       returns the diagonally trailing coefficient matrix
%
% The second output argument D in all cases returns the corresponding
% matrix or vector of trailing degrees. D is the same as the first output
% argument of the function TDEG.
%
% See also TDEG, POL/TCOEF, POL/TDEG, TSP/TCOEF, TSP/TDEG.

%          Author: J. Jezek  13-10-99
%          Copyright (c) 1999 by Polyx, Ltd.
%          $ Revision $  $ Date 15-Jun-2000 $
%
% Effect on other properties:
% T and D are standard Matlab matrices.

if nargin<1,
   error('Not enough input arguments.');
elseif nargin==1,
   string = 'mat';
elseif ~ischar(string),
   error('Invalid 2nd argument.');
end;

eval('[T,D] = tcoef(pol(A),string);','error(peel(lasterr));');

%end .. tcoef
