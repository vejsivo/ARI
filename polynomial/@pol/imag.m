function [Ai,Ar] = imag(A,tol);
%IMAG  Imaginary part of polynomial
%
% IMAG(A) is the imaginary part of A.
% [AI,AR] = IMAG(A) returns also the real part.
%
% See also POL/REAL.

%       Author(s): M. Hromcik, M. Sebek 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 06-Mar-1998 09:23:34   $
%       $Revision: 3.0 $  $Date: 28-Feb-2003  J.Jezek  tol  $

% Effect on other properties:
% Ar.u: UserData are deleted.

if nargin==2 & ~isempty(tol),
   if ~isa(tol,'double'),
      error('Invalid 2nd argument.');
   end;
end;

Ai = A;
Ai.c = imag(A.c);
Ai = pclear(Ai);

if nargout==2,
   Ar = A;
   Ar.c = real(A.c);
   Ar = pclear(Ar);
end;

%end .. @pol/imag
