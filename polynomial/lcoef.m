function [L,D] = lcoef(A,string)
%LCOEF  Leading coefficient of constant
%
% L = LCOEF(A)       default, the same as LCOEF(A,'mat')
% L = LCOEF(A,'mat') returns the leading coefficient matrix of A
% L = LCOEF(A,'ent') returns the scalar leading coefficients of
%                    the entries of A
% L = LCOEF(A,'row') returns the row leading coefficient matrix of A
% L = LCOEF(A,'col') returns the column leading coefficient matrix of A
% L = LCOEF(A,'dia') for symmetric matrix A,
%                     returns the diagonally leading coefficient matrix
%
% The second output argument D in all cases returns the corresponding
% matrix or vector of degrees. D is the same as the first output argument
% of the function DEG.
%
% See also DEG, POL/LCOEF, POL/DEG, TSP/LCOEF, TSP/DEG.

%       Author(s):  S. Pejchova, M. Sebek 17-3-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 22-Apr-1999 12:20:34   $
%       $Revision: 3.0 $  $Date: 15-Jun-2000 12:00:00   J. Jezek  $

% Effect on other properties:
% L and D are standard Matlab matrices.

if nargin<1,
   error('Not enough input arguments.');
elseif nargin==1,
   string='mat';
elseif ~ischar(string),   
  error('Invalid 2nd argument.');
end;

eval('[L,D] = lcoef(pol(A),string);','error(peel(lasterr));');

%end .. lcoef
