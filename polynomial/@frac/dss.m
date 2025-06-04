function sys = dss(F,T,tol)
%DSS   Convert fraction to Control Toolbox descriptor state space
%
% For proper fraction F, the command
%     SYS = DSS(F[,T[,TOL]])
% converts a system represented by the fraction to a Control
% System Toolbox LTI object in descriptor state space form.
%
% Improper fraction F cannot be converted to Control Toolbox
% descriptor state space form. Use ABCDE but matrix E will be
% singular, inacceptable for Control Toolbox.
%
% When T, the sampling time in the case of discrete-time
% systems, is not specified, it is taken from F.h . When
% F.h is empty, the default value 1 is taken. An optional
% tolerance may be specified in TOL. Its default value is
% the global zeroing tolerance.
% 
% See also RDF/ABCDE, LDF/ABCDE, MDF/ABCDE, SDF/ABCDE.

%     Author:  J. Jezek, 16-Feb-2000
%     Copyright(c) 2000 by Polyx, Ltd.
%     $ Revision $  $ Date 31-May-2000 $
%                   $ Date 25-Jul-2002 $
%                   $ Date 04-Oct-2009 $ M. Sebek: improper fractions
%                   allowed

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modyfied by M. Sebek 2009: improper fractions now allowed
%if issingular(E),
%   error('Fraction is improper, matrix E singular. Control Toolbox does not accept it.');
%end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Imp = 'Control Toolbox not found.';
switch F.v
case {'s','p'}
   eval('sys = dss(A,B,C,D,E);','disp(Imp)');
otherwise
   eval('sys = dss(A,B,C,D,E,T);','disp(Imp)');
end;

%end .. @frac/dss
