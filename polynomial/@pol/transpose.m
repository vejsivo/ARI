function At = transpose(A);
%TRANSPOSE (.')  Transpose polynomial
%                  A.' 
%
% AT = A.'  or  AT = TRANSPOSE(A) returns the non-conjugate transpose
% of the polynomial matrix A.
%
% See also POL/CTRANSPOSE.

%       Author(s): M. Hromcik, M. Sebek 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 06-Mar-1998 10:08:34   $

% Effect on other properties:
% At.u: UserData are deleted.
   
if nargin>1, error('Too many input arguments.'); end;
if nargout>1, error('Too many output arguments.'); end;

At.d = A.d;
At.s = fliplr(A.s);
At.c = permute(A.c,[2,1,3]);
At.v = A.v;
At.h = A.h;
At.u = [];
At.version = 3.0;
    
At = class(At,'pol');

%end .. @pol/transpose
 
