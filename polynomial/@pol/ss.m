function sys = ss(N,D,varargin)
%SS  Convert a left or right matrix fraction to an LTI object in state space form
%
% Given two polynomial matrices N and D with D nonsingular, the command
%     SYS = SS(N,D[,'l'|'r'][,T[,TOL]])
% converts a system represented in a left or right polynomial 
% matrix fraction form to a Control System Toolbox LTI object 
% in state space form.
%
% The 'l'-option is the default. The default value for T, the
% sampling time in the case of discrete-time systems, is either
% taken from N,D or it is 1. A tolerance value may be specified
% in TOL. Its default value is the global zeroing tolerance.
%
% See also: LTI2LMF, LTI2RMF, LMF2SS, SS2LMF.

%    Author: R.C.W. Strijbos, November 2, 1998.
%    Modified by D. Henrion, July 9, 1999.
%    Modified for 3.0 by J.Jezek, March 28, 2000.
%                        J.Jezek, July 11, 2000.
%    Copyright 1998-99 by Polyx, Ltd.


global PGLOBAL

% Default values
mfd = 'l';
st = [];
tol = PGLOBAL.ZEROING;
Var = PGLOBAL.VARIABLE;

eval('N = pol(N);', 'error(peel(lasterr));');
if nargin==1, D = pol(eye(N.s(1)));
else eval('D = pol(D);', 'error(peel(lasterr));');
end;

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

[tv,Var,N,D] = testvp(N,D);
if tv==2,
   if isempty(st),
      [tst,st,N,D] = testhp(N,D,Var);
      if ~tst,
         warning('Inconsistent sampling periods.');
      end;
   end;
   if strcmp(D.v,'z^-1'),
      dg = deg(D);
      N = shift(N,dg); D = rev(D,dg);
      D.v = 'z';
   else
      dg = deg(N);
      N = rev(N,dg); D = shift(D,dg);
      N.v = 'z';
   end;
   Var = 'z';
elseif ~tv,
   error('Inconsistent variables.');
end;

if isempty(st),
   [tst,st,N,D] = testhp(N,D,Var);
   if ~tst,
      warning('Inconsistent sampling periods.');
   end;
elseif length(st)~=1 | ~isreal(st) | st<0,
   error('Invalid sampling period.');
end;
if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

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
   a = 0; b = 0; c = 0; d = 0;  
   eval('[a,b,c,d] = lmf2ss(N,D,tol);', ...
      'error(peel(lasterr));');
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
   a = 0; b = 0; c = 0; d = 0;
   eval('[a,b,c,d] = rmf2ss(N,D,tol);', ...
      'error(peel(lasterr));');
end

if isa(d,'pol')
   error('Polynomials do not define proper transfer function.');
end

Imp='Control Toolbox not found.';
switch Var
case {'s','p'}
   eval('sys = ss(a,b,c,d);','disp(Imp)');
otherwise
   if isnan(st), st = [];
   end;
   eval('sys = ss(a,b,c,d,st);','disp(Imp)');
end;

%end .. @pol/ss
