function F = unsamp(G,varargin)
%UNSAMP     Unsample discrete time fraction
%
% Discrete time fraction G is supposed to be a result
% of sampling of continuous time fraction F with
% period T and phase TAU. The command
%     F  = UNSAMP(G,T,TAU)  restores the original
% fraction F.
%
% Fraction G, as a function of complex variable,
% must be proper, i.e. it must not have a pole for
% z=inf (z^-1=0, q=inf, d=0). Furthermore, it must
% have a zero point for z=0 (z^-1=inf, q=0, d=inf).
%
% The unsampling is not unique; the function F is
% selected so that imag parts of all poles are
% between -PI/T and PI/T. The "practically reasonable"
% case for control systems is: all poles of G are
% in the right halfplane of z. In such a case, the
% imag parts of poles F are between -PI/(2*T) and
% PI/(2*T).
%
% For real fraction G, the resulting fraction F is
% usually also real. However, when G has a real negative
% pole, F may happen to be complex.
%
% The arguments T,TAU >= 0 are optional, they may be
% omitted or given by []. Their default values are G.h
% and 0 . The variable of F ('s' or 'p') is given
% by optional argument VAR or, as a default, by
% PGLOBAL.CONTVAR.
%
% In case of zero sampling period T, only special G is
% allowed: it must be equal to  K*z/(z-1)  or  K*q/(q-1),
%  K/(1-z^-1),  K/(1-d) . The result is equal to  K/s
% or  K/p .
%
% An optional argument TOL may specify the zeroing
% tolerance used istead of the standard one. It must have
% a form of character string, e.g. '10^-6'. An optional
% argument MET may specify the method used. The methods
% are  'ss' - state space
%
% See also RDF/SAMP, LDF/SAMP, MDF/SAMP, SDF/SAMP.

%       Author: J.Jezek, 17-Oct-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 28-Dec-2000 $
%                     $ Date 25-Jul-2002 $
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

T = G.h; tau = 0; tol = PGLOBAL.ZEROING;
met = 'ss'; var = PGLOBAL.CONTVAR;
lv = length(varargin); argT = 1;
for i = 1:lv,
   arg = varargin{i};
   if isa(arg,'pol'),
      [vs1,vs2,vd] = size(arg);
      if all ([vs1,vs2,vd]==1) & all(arg.c(:,:)==[0 1]),
         arg = arg.v;
      else
         error(['Invalid ',nth(i+1),' argument.']);
      end;
   end;   
   if isa(arg,'char'),
      num = str2num(arg);
      if isa(num,'double') & ~isempty(num),
         if length(num)==1 & isreal(num) & num>=0 & num<=1,
            tol = num;
         else
            error(['Invalid ',nth(i+1),' argument.']);
         end;
      elseif ~isempty(arg),
         if strcmp(arg,'s') | strcmp(arg,'p'),
            var = arg;
         elseif strcmp(arg,'ss'),
            met = arg;
         else
            error(['Invalid ',nth(i+1),' argument.']);
         end;   
      end;
   elseif isa(arg,'double'),
      if ~isempty(arg),
         if argT,
            if length(arg)==1 & isreal(arg) & arg>=0,
               T = arg;
            else
               error('Invalid sampling period.');
            end;
         else
            if length(arg)==1 & isreal(arg) & arg>=0,
               tau = arg;
            else
               error('Invalid sampling phase.');
            end; 
         end;
      end;
      argT = 0;
   else
      error(['Invalid ',nth(i+2),' argument.']);
   end;
end;
      
if isempty(G) | all(all(G.num==0)),
   F = G; return;
end;

Gvar = G.v;
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
case '',
   error('Invalid 1st argument; has not required zero point.');
end;

savew = warning; warning off;
absterm = mvalue(G,0); warning(savew);
if any(any(absterm~=0)),
   error(['Invalid 1st argument; has not zero point ',text,'.']);
end;

if isempty(T),
   error('Sampling period missing.');
end;

if T==0,
   Gh = G.h;
   if isempty(Gh) | Gh~=0,
      num = shift(G.num,-1);
      if isa(num,'tsp'),
         error('Invalid 1st argument for zero old sampling period.');
      end;
      zz = z(num.h);
      num = num*(zz-1); den = G.den;
      [rG,cG] = size(G);
      Q = pol(ones(rG,cG)); R = Q;
      if isa(G,'rdf'),
         [Q,R] = rdiv(num,den);
      elseif isa(G,'ldf'),
         [Q,R] = ldiv(num,den);
      elseif isa(G,'mdf'),
         for i = 1:rG,
            for j = 1:cG,
               [Q(i,j),R(i,j)] = rdiv(num(i,j),den(i,j));
            end;
         end;
      elseif isa(G,'sdf'),
         for i = 1:rG,
            for j = 1:cG,
               [Q(i,j),R(i,j)] = rdiv(num(i,j),den);
            end;
         end;
      end;
      if deg(R)>=0 | deg(Q)>0,
         error('Invalid 1st argument for zero sampling period.');
      end;
   end;
   T = 1;
end;

if strcmp(met,'ss'),
   G.num = shift(G.num,-1);
   [AA,BB,CC,DD] = abcd(G,tol);
   A = logm(AA)/T; B = BB;
   if tau==0, C = CC;
   else C = CC*expm(-A*tau);
   end;
   
   strtol = num2str(tol);
   if isa(G,'rdf'),
      F = rdf(A,B,C,var,strtol);
   elseif isa(G,'ldf'),
      F = ldf(A,B,C,var,strtol);
   elseif isa(G,'mdf'),
      F = mdf(A,B,C,var,strtol);
   else
      F = sdf(A,B,C,var,strtol);
   end;
end;

%end .. @frac/unsamp

   
