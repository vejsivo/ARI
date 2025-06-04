function A = uplus(P);
%UPLUS (+)  Unary plus of polynomial
%            A = P, A = +P
%
% See also POL/UMINUS, POL/PLUS.

%       Author(s): M. Hromcik, M. Sebek 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 06-Mar-1998 10:13:34   $
%       $Revision: 3.0 $  $Date: 11-Aug-1999 12:00:00   $
%
% Effect on other properties: 
% A.u: UserData are deleted.
% A.version: set 3.0
 
A = P;
A.u = [];
A.version = 3.0;

%end .. @pol/uplus