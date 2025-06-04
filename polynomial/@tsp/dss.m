function sys = dss(F,T,tol)
%DSS   Convert two-sided polynomial to Control
%          Toolbox descriptor state space form
%
% For two-sided polynomial F, the command
%     SYS = DSS(F[T,[,TOL]])
% converts a system represented by F to a Control
% System Toolbox LTI object in descriptor state space form.
%
% If F contains positive powers of F then it cannot be
% converted to descriptor state space form. Use ABCDE
% but matrix E will bw singular, inacceptable for Control
% Toolbox.
%
% When T, the sampling period, is not specified, it is
% taken from F.h . Whwn F.h is empty, the default value 1
% is taken. An optional tolerance may be specified in TOL.
% Its default value is the global zeroing tolerance.
%
% See also TSP/ABCDE.

%      Author:  J. Jezek, 22-Feb-2003
%      Copyright(c) 2003 by Polyx, Ltd.

global PGLOBAL;

if nargin>=2,
   if ~isa(T,'double') | length(T)~=1 | ~isreal(T) | T<0,
      error('Invalid sampling period.');
   end;
   if nargin==3,
      if ~isa(tol,'double') | length(tol)~=1 | ~isreal(tol) | ...
            tol<0 | tol>1,
         error('Invalid tolerance.');
      end;
   else
      tol = PGLOBAL.ZEROING;
   end;
else
   if ~isempty(F.h), T = F.h;
   else T = 1;
   end;
   tol = PGLOBAL.ZEROING;
end;

[A,B,C,D,E] = abcde(F,tol);
if issingular(E),
   disp('Fraction is improper, matrix E singular.');
   error('Control Toolbox does not accept it.');
end;

Imp = 'Control Toolbox not found.';
switch F.v
case {'s','p'}
   eval('sys = dss(A,B,C,D,E);','disp(Imp)');
otherwise
   eval('sys = dss(A,B,C,D,E,T);','disp(Imp)');
end;

%end .. @tsp/dss
