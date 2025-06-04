function [D,L] = deg(A,string)
%DEG  Degree of constant
%
% D = DEG(A)       default, the same as DEG(A,'mat').
% D = DEG(A,'mat') returns the degree of matrix A.
% D = DEG(A,'ent') returns the matrix of degrees of the A entries.
% D = DEG(A,'row') returns the column vector of row degrees of A.
% D = DEG(A,'col') returns the row vector of column degrees of A.
% D = DEG(A,'dia') for symmetric matrix A,
%                   returns the vector of half diagonal degrees of A.
%
% The second output argument L in all cases returns corresponding
% matrix or vector of coefficients, the same as first output argument
% in function LCOEF.
%
% See also LCOEF, POL/DEG, POL/LCOEF, TSP/DEG, TSP/LCOEF.

%       Author(s):  S. Pejchova, M. Sebek 13-3-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 22-Apr-1999 12:18:34   $
%       $Revision: 3.0 $  $Date: 15-Jun-2000 12:00:00  J. Jezek  $

% Effect on other properties:
% D and L are standard Matlab matrices.

if nargin<1,
   error('Not enough input arguments.');
elseif nargin==1,
   string = 'mat';
elseif ~ischar(string),   
  error('Invalid command option.');
end;

eval('[D,L] = deg(pol(A),string);','error(peel(lasterr));');

%end .. deg
