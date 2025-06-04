function Y = value(F,X)
%VALUE      Value of two-sided polynomial
%
% The command  Y = VALUE(F,X)  returns matrix Y
% whose every entry Y(i,j) is the value of
% two-sided polynomial function F(i,j) of the
% argument X(i,j).
%
% Matrices F,X must have the same dimensions
% unless one of them is scalar. In such a case,
% the scalar operates with every entry of the
% other matrix.
%
% When F is empty or constant, the function does not
% depend on X. In such a case, X may be any scalar or
% any compatible matrix or it may be omitted.
%
% Matrix X or some entries of X may be infinite.
%
% See also POL/VALUE, TSP/MVALUE.

%     Author: J.Jezek, 09-Oct-2000
%     Copyright(c) 2000 by Polyx, Ltd.
%     $ Revision $  $ Date 31-Oct-2000 $
%                   $ Date 19-Dec-2000 $
%                   $ Date 27-Apr-2003 $

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

Fp = F.p; Fo = F.o;
Fnneg = nneg(F); Fnpos = npos(F);
for i = 1:sF1,
   for j = 1:sF2,
      Xij = X(i,j);
      if isinf(Xij),
         Y(i,j) = value(Fnneg(i,j),Xij);
      elseif Xij==0,
         Y(i,j) = value(Fnpos(i,j),inf);
      else
         Y(i,j) = Xij^Fo * value(Fp(i,j),Xij);
      end;
   end;
end;

%end .. tsp/value
