function N = neg(P)
%NEG  Negative part of a polynomial
%
% The command
%     N = NEG(P)
% for the polynomial P, returns the part with negative
% powers only. It is nonzero only if P is polynomial in z^-1.
%
% See also POL/POS, POL/NNEG, POL/NPOS.

%        Author:  J. Jezek  11-8-99
%        Copyright (c) 1999 by Polyx, Ltd. 
     
N = P;
if ~isempty(P.d),
   Pv = P.v;
   if isempty(Pv) | ~strcmp(Pv,'z^-1'),
      Ps = P.s; N = pol(zeros(Ps(1),Ps(2)));
   else
      U = P.c; U(:,:,1) = 0; N.c = U;
   end;
end;

%end .. @pol/neg
