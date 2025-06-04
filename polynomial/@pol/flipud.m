function Af = flipud(A);
%FLIPUD  Flip a polynomial matrix in up/down direction
%
% The command
%    FLIPUD(X) 
% returns X with its columns preserved and its rows flipped
% in the up/down direction.
%
% See also POL/FLIPLR, POL/ROT90.

%       Author(s): M. Hromcik, M. Sebek 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 06-Mar-1998 09:20:34   $
%       $Revision: 3.0 $  $Date: 24-Jun-2001  J.Jezek   $

% Effect on other properties:
% Af.u: UserData are deleted.
   
if nargin >1, error('Too many input arguments.'); end;
if nargout>1, error('Too many output arguments.'); end;

Af = A;

Af.c = flipdim(A.c, 1);
Af.u = [];
Af.version = 3.0;

%end .. @pol/flipud
