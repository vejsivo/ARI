function [Bl,Br] = polpart(B)
%POLPART   Polynomial matrix part extraction
%
% For square para-Hermitian two-sided polynomial matrix B,
% the command
%    BR = POLPART(B)
% extracts the half part of B such that
%    B = BR' + BR .
%
% For square two-sided polynomial matrix B, the command
%    [BL,BR] = POLPART(B)
% extracts left and right parts of B such that
%    B = BL' + BR .
%
% See also TSP/CTRANSP.

%     Author; J. Jezek, 25-Aug-2001
%     Copyright(c) 2001 by Polyx, Ltd.

[Bs1 Bs2] = size(B);
if Bs1 ~= Bs2,
   error('Matrix is not square.');
end;
if Bs1==0,
   Bl = B; Br = B; return;
end;

if nargout == 1,
   residue = norm(B-B');
   if residue > eps,
      warning(['POLPART: Two-sided polynomial matrix is not para-Hermitian. Residue = ' ...
            num2str(residue)]);
   end;
   Bl = nneg(B); Bl = Bl - Bl{0}*.5;
else
   Br = nneg(B); Br = Br - Br{0}*.5;
   Bl = npos(B); Bl = Bl - Bl{0}*.5;
   Bl = Bl';
end;

%end .. @tsp/polpart

   