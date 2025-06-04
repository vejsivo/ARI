function [Ar,Ai] = real(A,tol);
%REAL  Real part of polynomial
%
%  AR = REAL(A) is the real part of A.
%  [AR,AI] = REAL(A) returns also the imag part.
%
%  See also POL/IMAG.

%       Author(s): M. Hromcik, M. Sebek 27-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 06-Mar-1998 09:58:34   $
%       $Revision: 3.0 $  $Date: 28-Feb-2003  J.Jezek  tol  $

% Effect on other properties:
% Ar.u: UserData are deleted.

if nargin==2 & ~isempty(tol),
   if ~isa(tol,'double'),
      error('Invalid 2nd argument.');
   end;
end;

Ar = A;
Ar.c = real(A.c);
Ar = pclear(Ar);

if nargout== 2,
   Ai = A;
   Ai.c = imag(A.c);
   Ai = pclear(Ai);
end;

%end .. @pol/real
