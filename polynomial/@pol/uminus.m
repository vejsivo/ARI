function A = uminus(P);
%UMINUS (-)  Unary minus of polynomial
%            A = -P
%
% See also POL/UPLUS, POL/MINUS.

%       Author(s): M. Hromcik, M. Sebek 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 06-Mar-1998 10:12:34   $
%       $Revision: 3.0 $  $Date: 11-Aug-1999 12:00:00   $
%                         $Date: 06-Feb-2002            $

% Effect on other properties:
% A.version: set 3.0

A = P;
A.c = -P.c;
A.u = [];
A.version = 3.0;

%end .. @pol/uminus
