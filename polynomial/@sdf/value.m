function Y = value(F,X)
%VALUE    Value of scalar-denominator fraction
%
% For scalar-den fraction F, the command
%  Y = VALUE(F,X)  returns matrix Y whose
% every entry Y(i,j) is the value of function
% F(i,j) of the argument X(i,j). Matrices F,X
% must have the same dimensions unless one of
% them is scalar. In such a case, the scalar
% operates with every entry of the other matrix.
%
% When F is empty or constant, the function does not
% depend on X. In such a case, X may be any scalar or
% any compatible matrix or it may be omitted.
%
% Matrix X or some entries of X may be infinite.
%
% See also POL/VALUE, POL/MVALUE, SDF/MVALUE.

%       Author: J.Jezek, 17-Dec-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 14-Oct-2002 $
%                     $ Date 27-Apr-2003 $

ni = nargin;
if ni==2,
   if ~isa(X,'double') | ndims(X)~=2,
      error('Invalid 2nd argument.');
   end;
end;

sF = size(F);
if ni==2,
   sX = size(X);
   if any(sF~=sX),
      if all(sX==1),
         X = repmat(X,sF);
      elseif all(sF==1),
         F = repmat(F,sX); sF = sX;
      else
         error('Matrices not of the same dimensions.');
      end;
   end;
end;

sF1 = sF(1); sF2 = sF(2);
Y = zeros(sF1,sF2);
if sF1==0 | sF2==0, return;
end;
if ni==1,
   [FD,FC] = declass(F);
   if strcmp(FC,'double'),
      Y = FD; return;
   else
      error('Not enough input arguments.');
   end;
end;


N = F.frac.num; D = F.frac.den;
Nd = deg(N,'ent'); Dd = deg(D);
for i = 1:sF1,
   for j = 1:sF2,
      Xij = X(i,j);
      Ndij = Nd(i,j);
      if Ndij>=0,
         if isfinite(Xij),
            NN = N.c(i,j,Ndij+1);
            for k = Ndij:-1:1,
               NN = NN*Xij + N.c(i,j,k);
            end;
            DD = D.c(1,1,Dd+1);
            for k = Dd:-1:1,
               DD = DD*Xij + D.c(1,1,k);
            end;
         else
            NDdij = max(Ndij,Dd);
            if Ndij==NDdij,
               NN = N.c(i,j,NDdij+1);
            else
               NN = 0;
            end;
            if Dd==NDdij,
               DD = D.c(1,1,NDdij+1);
            else
               DD = 0;
            end;
         end;
         Y(i,j) = NN/DD;
      end;
   end;
end;

%end .. @sdf/value

