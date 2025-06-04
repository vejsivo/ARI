function [N,D] = lti2lmf(sys,tol)
%LTI2LMF  Convert LTI object to left polynomial matrix fraction.
%
% Given an LTI object SYS in Control System Toolbox format,
% the command
%    [N,D] = LTI2LMF(SYS[,TOL]])
% converts the LTI object to the left polynomial matrix fraction D^-1 *N.
%
% A tolerance value may be specified in TOL. Its default value 
% is the global zeroing tolerance.
%
% See also: SS, TF, ZPK, LTI2RMF.

%    Author: R.C.W. Strijbos, November 13, 1998.
%    Copyright 1998 by Polyx, Ltd.
%    $ Revision 3.0 $  $ Date 31-May-2000  J. Jezek $
%                      $ 28-Jul-2001 J.Jezek  sampl per $


global PGLOBAL;
eval('PGLOBAL.ZEROING;','painit;');

if nargin==0,
   error('Not enough input arguments.');
end;
if nargin == 2 & ~isempty(tol),
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end
else
   tol = PGLOBAL.ZEROING;
end

if ~isa(sys,'lti')
   error('Invalid 1st argument; not an LTI object.');
else
   if isa(sys,'tf')
      num = sys.num;
      den = sys.den;
      var = sys.v;
      per = sys.ts;
      [N,D] = tf2lmf(num,den,tol);
      if strcmp(var,'z^-1')
         [N,D] = reverse(N,D,'l',tol);
      end
      N.v = var; D.v = var;
      N.h = per; D.h = per;
   elseif isa(sys,'zpk')
      Z = sys.z;
      P = sys.p;
      K = sys.k;
      var = sys.v;
      per = sys.ts;
      [N,D] = zpk2lmf(Z,P,K,tol);
      if strcmp(var,'z^-1') | strcmp(var,'d')
	 [N,D] = reverse(N,D,'l',tol);
      end
      N.v = var; D.v = var;
      N.h = per; D.h = per;
   else
      [a1,b1,c1,d1,st] = ssdata(sys);
      [N,D] = ss2lmf(a1,b1,c1,d1,tol);
      if st==0
         N.v = 's';
         D.v = 's';
      else
         if st<0, st = [];
         end;
         N.v = 'z';
         D.v = 'z';
      end
      N.h = st; D.h = st;
   end
end

%end .. lti2lmf



