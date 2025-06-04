function [Z,K] = pol2root(P);
%POL2ROOT  Extract the zeros and gains of a polynomial matrix
%
% The command
%    [Z,K] = POL2ROOT(P) 
% returns a cell array Z of roots and an array K of gains of the
% entries of the polynomial matrix P.
%
% See also ROOT2POL, POL, ROOTS.

%	Author(s): M. Hromcik, M. Sebek 15-9-98
%	Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 1.0 $  $Date: 16-Sep-1998 10:28:34  $ 
%       $Revision: 2.0 $  $Date: 10-Oct-2000 10:30:00 - Version 3.0, M. Hromcik $

if nargin<1,
   error('Not enough input arguments.');
end;
eval('P = pol(P);', 'error(peel(lasterr));');
Ps = P.s;
if Ps(1)>0 & Ps(2)>0,
   Pc = flipdim(P.c, 3);
   Pd = P.d;
   if isinf(Pd), Pd = 0; Pc = zeros(Ps);
   end;
   Pdeg = deg(P,'ent');
   Pdeg(isinf(Pdeg)) = 0;
   Pdeg(isnan(Pdeg)) = 0;
end;

Z = cell(Ps);
K = zeros(Ps);

if Ps(1)>0 & Ps(2)>0,
   for i = 1:Ps(1),
      for j = 1:Ps(2),
         Pij = Pc(i,j,:);
         K(i,j) = Pij(Pd-Pdeg(i,j)+1);
         Z{i,j} = roots(Pij);
      end; 
   end;   
end;

%end .. pol2root
