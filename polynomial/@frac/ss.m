function sys = ss(F,T,tol)
%SS    Convert fraction to Control Toolbox state space
%
% For proper fraction F, the command
%     SYS = SS(F[,T[,TOL]])
% converts a sytem represented by the fraction to a Control
% System Toolbox LTI object in state space form.
%
% When T, the sampling time in the case of discrete-time
% systems, is not specified, it is taken from F.h . When
% F.h is empty, the default value 1 is taken. An optional
% tolerance may be specified in TOL. Its default value is
% the global zeroing tolerance.
% 
% See also RDF/ABCD, LDF/ABCD, MDF/ABCD, SDF/ABCD.

%     Author:  J. Jezek, 16-Feb-2000
%     Copyright(c) 2000 by Polyx, Ltd.
%     $ Revision $  $ Date 31-May-2000 $
%                   $ Date 25-Jul-2002 $

global PGLOBAL;

Fcl = class(F);
if strcmp(Fcl,'frac'),
   error('Invalid 1st argument.');
end;

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
   if ~isempty(F.h), T = F.h;
   else T = 1;
   end;
   tol = PGLOBAL.ZEROING;
end;

[A,B,C,D] = abcd(F,tol);
if isa(D,'pol')
   error('Fraction is not proper.');
end;

Imp = 'Control Toolbox not found.';
switch F.v
case {'s','p'}
   eval('sys = ss(A,B,C,D);','disp(Imp)');
otherwise
   eval('sys = ss(A,B,C,D,T);','disp(Imp)');
end;

%end .. @frac/ss
