function sys = tf(N,D,varargin)
%TF  Convert a left or right matrix fraction to an LTI object in transfer function form.
%
% Given two polynomial matrices N and D, the command
%    SYS = TF(N,D)   or   SYS = TF(N)   or
%    SYS = TF(N,D [,'l' | 'r'] [,T [,TOL]])
% converts a system represented in a left or right polynomial matrix 
% fraction form to a Control System Toolbox LTI object in transfer 
% function form.
%
% If D is not present, a unit matrix of corresponding size is assumed.
% The 'l'-option is the default. The default value for T, the sampling 
% time in the case of discrete-time systems, is either taken from N,D 
% ot it is 1. An optional tolerance may be specified in TOL. Its
% default value is the global zeroing tolerance.
%
% See also: LTI2LMF, LTI2RMF, LMF2TF, TF2LMF.

%    Author: R.C.W. Strijbos, November 13, 1998.
%    Copyright 1998 by Polyx, Ltd.
%    Modified by Rens Strijbos on May 10, 2000
%         $ Revision 3.0 $  $ Date 28-Mar-2000  J.Jezek  $
%                           $ Date 11-Jul-2000  J.Jezek  $ 

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

if strcmp(Var, 'z^-1') | strcmp(Var, 'd'),
   [N,D] = reverse(N,D,'r',tol);
   N.v = 'z'; D.v = 'z'; Var = 'z';
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

if isempty(N) | isempty(D),
   sys = tf(zpk(N,D));
   return;
end;

switch mfd
case 'l'
   [num,den] = lmf2tf(N,D,tol);
case 'r'
   [num,den] = rmf2tf(N,D,tol);
end

Imp = 'Control Toolbox not found.';
switch Var
case {'s','p'}
   eval('sys = tf(num,den);','disp(Imp)');
case {'z','q'}
   if isnan(st), st = [];
   end;
   eval('sys = tf(num,den,st,''variable'',''z'');','disp(Imp)');
otherwise
   if iscell(num)
      [nr,nc]=size(num);
      for i = 1:nr
         for j = 1:nc
            num{i,j} = fliplr(num{i,j});
         end
      end
   else
      num = fliplr(num);
   end
   if iscell(den)
      [nr,nc]=size(den);
      for i = 1:nr
         for j = 1:nc
           den{i,j} = fliplr(den{i,j});
         end
      end
   else
      den = fliplr(den);
   end
   
   eval('sys = tf(num,den,st,''variable'',''z^-1'');','disp(Imp)');
end;

%end .. @pol/tf

