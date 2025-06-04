function F = unsamph(G,varargin)
%UNSAMPH     Un- sample and hold discrete time fraction
%
% Discrete time fraction G is supposed to be a result
% of sampling and holding of continuous time fraction
% F with period T and phase TAU. The command
%     F = UNSAMPH(G,T,TAU)
% restores the original fraction F.
%
% Fraction G, as a function of complex variable,
% must be proper, i.e. it must not have a pole for
% z=inf (z^-1=0, q=inf, d=0). Furthermore, it must
% not have a pole for z=0 (z^-1=inf, q=0, d=inf).
%
% The arguments T,TAU are optional, they may be
% omitted or given by []. The defaults are G.h and 0 .
% The variable of F ('s' or 'p') is given by optional
% argument VAR or, as a default, by PGLOBAL.CONTVAR.
%
% An optional argument TOL may specify the zeroing
% tolerance used istead the standard one. It must have
% a form of character string, e.g. '10^-6'. An optional
% argument MET may specify the method used. The methods
% are  'ss' - state space
%
% See also FRAC/SAMPH, FRAC/UNSAMP.

%       Author: J.Jezek, 17-Oct-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 25-Dec-2000 $
%                     $ Date 25-Jul-2002 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $
%                     $ Date 28-Feb-2003  warning  $

global PGLOBAL;

Gcl = class(G);
if strcmp(Gcl,'frac'),
   error('Invalid 1st argument.');
end;
if ~isa(G,'frac'),
   error('Some argument but not 1st is invalidly fraction.');
end;
T = G.h;

tol = PGLOBAL.ZEROING;
lv = length(varargin);
if lv>=1,
   arg = varargin{1};
   if ~isempty(arg),
      if ~isa(arg,'double') | length(arg)~=1 | ...
            ~isreal(arg) | arg<=0,
         error('Invalid sampling period.');
      end;
      T = arg;
   end;

   for i = 2:lv,
      arg = varargin{i};
      if isa(arg,'char'),
         num = str2num(arg);
         if isa(num,'double') & ~isempty(num),
            if length(num)==1 & isreal(num) & num>=0 & num<=1,
               tol = num;
            else
               error(['Invalid ',nth(i+1),' argument.']);
            end;
         end;
      end;
   end;
end;

Gvar = G.v;
if isempty(Gvar) | all(all(G.num==0)),
   F = G; return;
end;
if strcmp(Gvar,'s') | strcmp(Gvar,'p'),
   error('Invalid 1st argument; must be discrete time.');
end;
   
if ~isproper(G),
   error('Invalid 1st argument; is not proper.');
end;   
switch Gvar,
case 'z', text = 'z=0';
case 'q', text = 'q=0'; G.v = 'z';
case 'z^-1', text = 'z^-1=inf'; G = reverse(G);
case 'd', text = 'd=inf'; G.v = 'z^-1'; G = reverse(G);
end;

if T==0,
   error('Sampling period is zero.');
end;
zz = z(T);
G.h = T; G.num.h = T; G.den.h = T;
savew = warning; warning off;
absterm = mvalue(G,0); warning(savew);
if ~all(all(isfinite(absterm))),
   error(['Invalid 1st argument; has pole ',text,'.']);
end;

G.num = G.num*zz; G.den = G.den*(zz-1);
eval('F = unsamp(G,T,varargin{2:lv});', ...
   'error(peel(lasterr));');
var = pol([0 1],1,F.v);
F.num = F.num*var;

F.c = 'cop?'; F.r = 'red?'; F.p = 'prop';
F.tc = []; F.tp = [];
if strcmp(PGLOBAL.COPRIME,'cop'),
   F = coprime(F,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   F = reduce(F,tol);
else
   F = smreduce(F);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'),
   F = defract(F);
end;

%end .. @frac/unsamph
