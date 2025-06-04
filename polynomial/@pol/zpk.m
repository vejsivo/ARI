function sys = zpk(N,D,varargin)
%ZPK  Convert a left or right matrix fraction to an LTI object in zero-pole-gain form
%
% Given two polynomial matrices N and D, the command
%    SYS = ZPK(N,D)   or   SYS = ZPK(N)   or
%    SYS = ZPK(N,D[,'l'|'r'][,T[,TOL]])
% converts a system represented in left or right polynomial 
% matrix fraction form to a Control System Toolbox LTI object 
% in zero-pole-gain form.
%
% If D is not present, a unit matrix of corresponding size is assumed.
% The 'l'-option is the default and the default value for T, 
% the sampling time in the case of discrete-time systems, is
% either taken from N,D or it is 1. 
% An optional tolerance may be specified in TOL. Its default 
% value is the global zeroing tolerance.
%
% See also: LTI2LMF, LTI2RMF, LMF2ZPK, ZPK2LMF.

%    Author: R.C.W. Strijbos, November 13, 1998.
%    Copyright 1998 by Polyx, Ltd.
%        $ Revision 3.0 $  $ Date: 28-Mar-2000  J.Jezek $
%    Modified by Rens Strijbos on May 10, 2000
%        $ Revision $      $ Date: 11-Jul-2000  J.Jezek $

global PGLOBAL

% Default values
mfd = 'l';
st = 1;
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
   if strcmp(Var,'z^-1') | strcmp(Var,'d')
      [N,D] = reverse(N,D,'l',tol);
      N.v = 'z'; D.v = 'z'; Var = 'z';
   end
   [Z,P,K] = lmf2zpk(N,D,tol);
case 'r'
   if strcmp(Var,'z^-1') | strcmp(Var,'d')
      [N,D] = reverse(N,D,'r',tol);
      N.v = 'z'; D.v = 'z'; Var = 'z';
   end   
   [Z,P,K] = rmf2zpk(N,D,tol);
end

Imp = 'Control Toolbox not found.';
switch Var
case {'s','p'}
   eval('sys = zpk(Z,P,K);','disp(Imp)');
otherwise
   eval('sys = zpk(Z,P,K,st);','disp(Imp)');
end;

%end .. @pol/zpk
