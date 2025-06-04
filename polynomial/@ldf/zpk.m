function sys = zpk(F,T,tol)
%ZPK  Convert left-den fraction to LTI object
%     in zero-pole-gain form
%
% For left-den fraction F, the command
%    SYS = TF(F[,T,[TOL]])
% converts a system represented by the fraction to a Control
% System Toolbox object in zero-pole-gain form.
%
% When T, the sampling time in the case of discrete-time
% systems, is not specified, it is taken from F.h . When
% F.h is empty, the default value 1 is taken. An optional
% tolerance may be specified in TOL. Its default value is
% the global zeroing tolerance.
%
% For backward conversion, see LDF/LDF.

%     Author:  J. Jezek, 16-Feb-2000
%     Copyright(c) 2000 by Polyx, Ltd.
%     $ Revision $  $ Date 31-May-2000 $
%                   $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin>=2,
   if ~isa(T,'double') | length(T)~=1 | ~isreal(T) | T<0,
      error('Invalid sampling period.');
   end;
   if nargin==3,
      if ~isa(tol,'double') | length(tol)~=1 | ...
            ~isreal(tol) | tol<0 | tol>1,
         error('Invalid tolerance.');
      end;
   else
      tol = PGLOBAL.ZEROING;
   end;
else
   if ~isempty(F.frac.h), T = F.frac.h;
   else T = 1;
   end;
   tol = PGLOBAL.ZEROING;
end;

[Z,P,K] = lmf2zpk(F.frac.num,F.frac.den,tol);

Imp = 'Control Toolbox not found.';
switch F.frac.v
case {'s','p'}
   eval('sys = zpk(Z,P,K);','disp(Imp)');
otherwise
   eval('sys = zpk(Z,P,K,T);','disp(Imp)');
end;

%end .. @ldf/zpk


