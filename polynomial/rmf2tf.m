function [num,den] = rmf2tf(N,D,tol)
%RMF2TF  Convert a right matrix fraction to a Control System Toolbox transfer function
%
% Given two polynomial matrices N and D, where D is square and
% nonsingular with the same number of rows as N, the function
%    [NUM,DEN] = RMF2TF(N, D [,TOL])
% returns the transfer function H = N*D^{-1} in Control System 
% Toolbox format.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: LMF2TF, TF2RMF.

%    Authors: D. Henrion, M. Sebek, R.C.W. Strijbos, November 13, 1998.
%    Updated to 3.0 by D. Henrion, July 25, 2000.
%                   by J. Jezek, Aug 20, 2001,  arg checking
%    Updated by D. Henrion, November 25, 2002, scalar case
%    Copyright 1998-2002 by Polyx, Ltd.

%    Numerator and denominator entries in H are computed from
%    the adjoint matrix and determinant of D, i.e. H = N*adj(D)/det(D).
%    Common factors are cancelled with macros GRD and XAB.

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');
  
switch nargin
case 2
   tol = PGLOBAL.ZEROING;
case 3
   if isempty(tol),
      tol = PGLOBAL.ZEROING;
   elseif ~isa(tol,'double') | length(tol) > 1 ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.')
   end
otherwise
   error('Not enough input arguments.');
end

eval('N = pol(N); D = pol(D);', ...
   'error(peel(lasterr));');
eval('[N,D] = testdnd(N,D);', ...
   'error(peel(lasterr));');
Var = ''; h = 0;
eval('[Var,h,N,D] = testvhnd(N,D);', ...
   'error(peel(lasterr));');
   
if strcmp(Var, 'z^-1') | strcmp(Var, 'd'),
   [N,D] = reverse(N,D,'r',tol);
   N.v = 'z'; D.v = 'z';
end;

[rN,cN] = size(N);
[rD,cD] = size(D);

if max([rN cN rD cD]) > 1,
  
 % matrix transfer function
 [AdjD, detD] = adj(D, tol);
 if pzer(detD) == 0
   error('Denominator matrix is singular.');
 end

 Num = N * AdjD;
 num = pol(zeros(rN,cN)); den = num;
 for i = 1:rN,
   for j = 1:cN,

      numerator = Num(i,j);
      factor = grd(numerator, detD, tol);
      polN = xab(factor, numerator, tol);
      polD = xab(factor, detD, tol);

      if isnan(polN) | isnan(polD),
         error('Factorization failed. Try to modify tolerance.');
      end;
      num(i,j) = polN;
      den(i,j) = polD;
   end;
 end;

else

 % scalar transfer function
 num = N;
 den = D;

end;

num = pzer(num);
den = pzer(den);

num = pol2mat(num);
den = pol2mat(den);

%end .. rmf2tf

