function Y = mvalue(F,X)
%MVALUE      Matrix value of polynomial
%
% The command  Y = MVALUE(F,X)  returns value F(X)
% of polynomial function F of the argument X.
%
% If F is scalar polynomial then X may be square
% matrix; otherwise it must be scalar.
%
% When F is empty or constant, the function does not
% depend on X. In such a case, X may be any scalar or
% it may be omitted.
%
% Matrix X or some entries of X may be infinite.
%
% See also POL/VALUE.

%       Author: J.Jezek, 05-Oct-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 16-Dec-2000 $
%                     $ Date 27-Apr-2003 $

ni = nargin;
if ni==2,
   if ~isa(X,'double') | ndims(X)~=2,
      error('Invalid 2nd argument.');
   end;
end;

if length(F)==1,
   if ni==2,
      sX = size(X);
      if sX(1)~=sX(2),
         error('Invalid 2nd argument; must be square matrix.');
      end;
   end;
   Fd = F.d;
   if Fd<0,
      if ni==2, Y = zeros(sX);
      else Y = 0;
      end;
      return;
   end;
   Y = F.c(Fd+1);
   if ni==1,
      if Fd==0, return;
      else error('Not enough input arguments.');
      end;
   end;
   I = eye(sX);
   Y = Y*I;
   if Fd>=1,
     for k = Fd:-1:1,
         Y = Y*X + F.c(k)*I;
     end;
   end;
  
else
   if ni==2,
      if length(X)~=1,
         error('Invalid 2nd argument; must be scalar.');
      end;
   end;
   sF = size(F); sF1 = sF(1); sF2 = sF(2);
   Fd = F.d;
   if sF1==0 | sF2==0 | Fd<0,
      Y = zeros(sF); return;
   end;
   if ni==1,
      if Fd==0, Y = F.c(:,:,1); return;
      else error('Not enough input arguments.');
      end;
   end;
   Y = F.c(:,:,Fd+1);
   if Fd>=1,
      for k = Fd:-1:1,
         Y = Y*X + F.c(:,:,k);
      end;
   end;
   
end;

%end .. @pol/mvalue
