function Af = fliplr(A);
%FLIPLR  Flip a polynomial matrix in left/right direction
%
% The command
%    FLIPLR(A) 
% returns A with its rows preserved and its columns flipped in the 
% left/right direction.              
%
% See also POL/FLIPUD, POL/ROT90.

%       Author(s): M. Hromcik, M. Sebek 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 06-Mar-1998 09:18:34   $
%       $Revision: 3.0 $  $Date: 24-Jun-2001  J.Jezek   $

% Effect on other properties:
% Af.u: UserData are deleted.
   
if nargin >1, error('Too many input arguments.'); end;
if nargout>1, error('Too many output arguments.'); end;

Af = A;

Af.c = flipdim(A.c, 2);
Af.u = [];
Af.version = 3.0;

%end .. @pol/fliplr
