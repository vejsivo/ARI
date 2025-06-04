function Acj = conj(A);
%CONJ   Complex conjugate of polynomial
%
% The commmand
%    CONJ(A) 
% returns the complex conjugate of A, that is
%    CONJ(A) = REAL(A) - i*IMAG(A)
%
% See also POL/CTRANSPOSE, POL/REAL, POL/IMAG.

%  	Author(s): M. Hromcik, M. Sebek 16-2-98
%	   Copyright (c) 1998 by Polyx, Ltd.
%     $Revision: 2.0 $  $Date: 09-Mar-1998 13:48:34   $

% Effect on other properties:
% Acj.u: UserData are deleted.
   
if nargin > 1,
   error('Too many input arguments.'); 
elseif nargout > 1,
   error('Too many output arguments.');
end; 
 
Acj.d = A.d;
Acj.s = A.s;
Acj.c = conj(A.c);
Acj.v = A.v;
Acj.h = A.h;
Acj.u = [];
Acj.version = 3.0;

Acj = class(Acj, 'pol');

%end .. @pol/conj
