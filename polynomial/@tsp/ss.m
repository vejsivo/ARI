function sys = ss(F,T,tol)
%SS    Convert two-sided polynomial to Control
%          Toolbox state space form
%
% For two-sided polynomial F, the command
%     SYS = SS(F[T,[,TOL]])
% converts a system represented by F to a Control
% System Toolbox LTI object in state space form.
%
% If F contains positive powers of F then it does not
% represent a proper transfer function and cannot be
% converted to Control Toolbox state space form. Use ABCD
% but matrix D will bw polynomial, inacceptable for Control
% Toolbox.
%
% When T, the sampling period, is not specified, it is
% taken from F.h . Whwn F.h is empty, the default value 1
% is taken. An optional tolerance may be specified in TOL.
% Its default value is the global zeroing tolerance.
%
% See also TSP/ABCD.

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

[A,B,C,D] = abcd(F,tol);
if ~isa(D,'double'),
   error('Polynomials do not define proper transfer function.');
end;

Imp = 'Control Toolbox not found.';
switch F.v
case {'s','p'}
   eval('sys = ss(A,B,C,D);','disp(Imp)');
otherwise
   eval('sys = ss(A,B,C,D,T);','disp(Imp)');
end;

%end .. @tsp/ss


