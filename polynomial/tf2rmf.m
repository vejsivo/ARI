function [N,D] = tf2rmf(num, den, tol)
%TF2RMF  Convert a Control System Toolbox transfer function to a right matrix fraction
%
% Given two cell arrays NUM and DEN of the same dimensions, the function
%    [N,D] = TF2RMF(NUM, DEN [,TOL])
% converts the rational transfer matrix (in MATLAB format) defined by
% NUM and DEN to the right coprime polynomial matrix fraction N*D^-1.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
%    See also: TF2LMF, RMF2TF.

%    Authors: D. Henrion , M. Sebek, R.C.W. Strijbos, November 13, 1998.
%    Copyright 1998 by Polyx, Ltd.
%    Modified by J.Jezek 24-Jul-2001, empty case

% The macro calls LRM to extract the LCM of the denominators of the TF
% and then calls LMF2RMF to retrieve a coprime RMF.


global PGLOBAL;
eval('PGLOBAL.ZEROING;','painit;');

switch nargin
case {0,1}
   error('Not enough input arguments.');
case 2
   tol = PGLOBAL.ZEROING;
case 3
   if isempty(tol),
      tol = PGLOBAL.ZEROING;
   elseif ~isa(tol,'double') | length(tol) > 1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.')
   end
end
 
eval('num = mat2pol(num); den = mat2pol(den);', ...
   'error(peel(lasterr));');

[rN,cN] = size(num); [rD,cD] = size(den);
if rN~=rD | cN~=cD
   error('Matrices not of the same dimensions.');
end

% Build polynomial matrix of numerators

LCMT = pol(zeros(1,rN));
if rN>0,
   for i = 1:rN
      if cN==0,
         LCM = 1;
      else
         dencell = cell(1,cN);
         for j = 1:cN
            dencell{j} = den(i,j); %Extract least common multiple of denominators
         end
         LCM = lrm(dencell{:}, tol);
         if LCM == 0
            error('Invalid denominator.');
         end
      end;
      LCMT(i) = LCM;
      if cN>0,
         for j = 1:cN,
            num(i,j) = num(i,j) * rdiv(LCM, den(i,j), tol);
         end;
      end;   
   end;
end;

[N,D] = lmf2rmf(num, diag(LCMT), tol);

%end .. tf2rmf
