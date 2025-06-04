function NN = nneg(P)
%NNEG  Nonnegative part of a polynomial

% The command
%     NN = NNEG(P)
% for the polynomial P, returns the part with nonnegative
% powers only. If P is the polynomial in z^-1, the result
% is only the absolute term, if any.
%
% See also POL/NPOS, POL/NEG, POL/POS.

%        Author:  J. Jezek  11-8-99
%        Copyright (c) 1999 by Polyx, Ltd.

NN = P; 
if ~isempty(P.v) & strcmp(P.v,'z^-1'),
   U = P.c; Ps = P.s;
   V = U(Ps(1),Ps(2),1); NN = pol(V);
end;

%end .. @pol/nneg
