function Q = coprime(R,tol)
%COPRIME    Make matrix-den fraction coprime
%
% For matrix-den fraction R, the command  Q = COPRIME(R)  returns
% matrix-den fraction Q, which is equal to R but modified so that
% all its elements are coprime fractions.
%
% An optional input argument TOL may specify a zeroing tolerance 
% to be used instead the standard one.
%
% See also MDF/REDUCE.

%       Author:  J. Jezek  06-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 26-Apr-2000 $
%                     $ Date 06-Feb-2001 $
%                     $ Date 30-Jul-2001  bug 'red?' $
%                     $ Date 06-Oct-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin==1 | isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double') | length(tol)~=1 | ...
      ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

Q = R;
if strcmp(R.frac.c,'cop') & R.frac.tc==tol,    % quick exit
   return;
end;

Rs1 = R.frac.s(1); Rs2 = R.frac.s(2);
if Rs1>0 & Rs2>0,
   for i = 1:Rs1,
      for j = 1:Rs2,
         [Q.frac.num(i,j),Q.frac.den(i,j)] = ...
            axby0(-R.frac.den(i,j),R.frac.num(i,j),tol);
      end;
   end;
   Q.frac.h = R.frac.h;
   props(Q,'red?');
end;
props(Q,'cop',tol,R.frac.p,R.frac.tp);

%end .. @mdf/coprime
