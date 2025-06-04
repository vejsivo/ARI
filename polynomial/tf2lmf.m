function [N,D] = tf2lmf(num, den, tol)
%TF2LMF  Convert a Control System Toolbox transfer function to a left matrix fraction
%
% Given two cell arrays NUM and DEN of the same dimensions, the function
%    [N,D] = TF2LMF(NUM, DEN [,TOL])
% converts the rational transfer matrix (in MATLAB format) defined by
% NUM and DEN to the left coprime polynomial matrix fraction D^-1 * N.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: TF2RMF, LMF2TF.

%    Authors: D. Henrion , M. Sebek, R.C.W. Strijbos, November 13, 1998.
%    Copyright 1998 by Polyx, Ltd.
%    Modified by J.Jezek 24-Jul-2001, empty case

% The macro calls LRM to extract the LCM of the denominators of the TF
% and then calls RMF2LMF to retrieve a coprime LMF

global PGLOBAL;
eval('PGLOBAL.ZEROING;','painit');

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

LCMT = pol(zeros(1,cN));
if cN>0,
   for i = 1:cN
      if rN==0,
         LCM = 1;
      else
         dencell = cell(rN,1);
         for j = 1:rN
            dencell{j} = den(j,i); %Extract least common multiple of denominators
         end
         LCM = lrm(dencell{:}, tol);
         if LCM == 0
            error('Invalid denominator.');
         end
      end;
      LCMT(i) = LCM;
      if rN>0,
         for j = 1:rN,      
            num(j,i) = num(j,i) * rdiv(LCM, den(j,i), tol);
         end;
      end;
   end;
end;

[N,D] = rmf2lmf(num, diag(LCMT), tol);

%end .. tf2lmf

