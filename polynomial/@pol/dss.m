function sys = dss(N,D,varargin)
%DSS  Convert a left or right matrix fraction to an LTI object
% in descriptor state space form
%
% Given two polynomial matrices N and D with D nonsingular,
% the command   SYS = DSS(N,D)   or   SYS = DSS(N)   or
%     SYS = DSS(N,D[,'l'|'r'][,T[,TOL]])
% converts a system represented in a left or right polynomial 
% matrix fraction form to a Control System Toolbox LTI object 
% in state space form.
%
% IF D is not present, a unit matrix of corresponding size is assumed.
% The 'l'-option is the default. The default value for T, the
% sampling time in the case of discrete-time systems, is 1. 
% A tolerance value may be specified in TOL. Its default value 
% is the global zeroing tolerance.
%
% See also: LTI2LMF, LTI2RMF, LMF2SS, SS2LMF.

%    Author: J.Jezek, March 28, 2000.
%    Copyright(c) 2000 by Polyx, Ltd.
%    $ Revision $  $ Date 28-Feb-2003 $

global PGLOBAL

% Default values
mfd = 'l';
st = 1;
tol = PGLOBAL.ZEROING;
Var = PGLOBAL.VARIABLE;

eval('N = pol(N);', 'error(peel(lasterr));');
if nargin==1,
   D = pol(eye(N.s(1)));
else
   eval('D = pol(D);', 'error(peel(lasterr));');
end;

if isempty(N.v) & ~isempty(D.v)
   Var = D.v;
end
if isempty(D.v) & ~isempty(N.v)
   Var = N.v;
end
if ~isempty(N.v)
   Var = N.v;
end

if ~strcmp(N.v,D.v) & ~(isempty(N.v) | isempty(D.v))
   error('Inconsistent variables.');
end

if nargin >= 3
   a=varargin{1};
   if isa(a,'double') & length(a)==1
      st = a;
   elseif isa(a,'char') & (strcmp(a,'l') | strcmp(a,'r'))
      mfd = a;
   else
      error('Invalid 3rd argument.');
   end
   if nargin >= 4
      b=varargin{2};
      if ~(isa(b,'double') & length(b)==1)
         error('Invalid 4th argument.');
      end
      if isa(a,'double')
         tol = b;
      else
         st = b;
      end
      if nargin == 5
         c = varargin{3};
	 if ~(isa(c,'double') & length(c)==1)
	    error('Invalid 5th argument.');
	 else
	    tol = c;
         end
      end
   end
end

switch mfd
case 'l'
   % Check if D is row reduced
   [rD, cD] = size(D);
   rowD = lcoef(D, 'row');
   S = svd(rowD);
   if min(S) < rD * max(S) * tol,
     % Perform automatic row reduction
     [D,rkD,U] = rowred(D,tol);
     N = U*N;
   end;
   [a,b,c,d,e] = lmf2dss(N,D,tol);
case 'r'
   % Check if D is column reduced
   [rD, cD] = size(D);
   colD = lcoef(D, 'col');
   S = svd(colD);
   if min(S) < cD * max(S) * tol,
     % Perform automatic column reduction
     [D,rk,U] = colred(D,tol);
     N = N*U;
   end;
   [a,b,c,d,e] = rmf2dss(N,D,tol);
end

if issingular(e),
   disp('Fraction is improper, matrix E singular.');
   error('Control Toolbox does not accept it.');
end;

Imp = 'Control Toolbox not found';
switch Var
case {'s','p'}
 eval('sys = dss(a,b,c,d,e);','disp(Imp)');
otherwise
 eval('sys = dss(a,b,c,d,e,st);','disp(Imp)');
end;

%end .. @pol/dss
