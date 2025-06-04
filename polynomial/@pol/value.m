function Y = value(F,X)
%VALUE       Value of polynomial
%
% The command  Y = VALUE(F,X)  returns matrix Y
% whose every entry Y(i,j) is the value of
% polynomial function F(i,j) of the argument X(i,j).
% 
% Matrices F,X must have the same dimension unless
% one of them is scalar. In such a case, the scalar
% operates with every entry of the other matrix.
%
% When F is empty or constant, the function does not
% depend on X. In such a case, X may be any scalar or
% any compatible matrix or it may be omitted.
%
% Matrix X or some entries of X may be infinite.
%
% See also POL/MVALUE.

%       Author: J.Jezek, 09-Oct-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 16-Dec-2000 $
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
         F = repmat(F,sX); sF = size(F);
      else
         error('Matrices not of the same dimensions.');
      end;      
   end;
end;

Y = zeros(sF);
Fd = deg(F);
if  any(sF==0) | Fd<0, return;
end;

Y = F.c(:,:,Fd+1);
if Fd>=1,
   if ni<2,
      error('Not enough input arguments.');
   end;
   for k = Fd:-1:1,
      Y = Y.*X + F.c(:,:,k);
   end;
end;

%end .. @pol/value
