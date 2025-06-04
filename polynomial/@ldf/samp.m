function G = samp(F,varargin)
%SAMP   Sample continuous time left-den fraction
%
% Let F be a left-denominator fraction in variable 's' or 'p'.
% The command  G = SAMP(F,T,TAU) where scalars T,TAU are
% sampling period, T>=0, and sampling phase,  TAU >=0 ,
% returns proper left-denominator fraction G in variable 'z','q',
% 'z^-1' or 'd', whose Laurent series equals to the inverse
% Laplace transform of F, sampled with period T and phase TAU,
% i.e. in sampling points  k*T + TAU, k = 0,1,2...
%
% For example:
%    F = (s+a)\1
%    inverse Laplace = exp(-a*t)
%    
%    SAMP(F,T,0) = (z-exp(-a*T))\z
%    laurent = 1 + exp(-a*T)*z^-1 + exp(-a*2T)*z^-2 + ...
%
%    SAMP(F,T,TAU) = (z-exp(-a*T))\exp(-a*TAU)*z
%    laurent = exp(-a*TAU) + exp(-a*(T+TAU))*z^-1 
%                          + exp(-a*(2T+TAU))*z^-2 + ...
%
% TAU can also be, instead of scalar, a vector. In this case,
% the results form a cell vector. All cells have the
% same denominator.
%
% The arguments T,TAU are optional, the defaults being
%  T = 1,  TAU = 0 .
% When TAU is zero, the fraction F must be strictly proper
% (otherwise, in time domain, sampling of a Dirac impulse
% would be attempted).
%
% The variable of G is given by optional argument VAR or, as a
% default, by PGLOBAL.DISCRVAR .
%
% An optional argument TOL may specify a zeroing tolerance to be
% used instead of the standard one. It must have a form of character
% string, e.g. '10^-6' .
%
% An optional argument MET may specify the method to be used.
% The methods are 'ss' - state space
%
% See also LDF/LAURENT, FRAC/LAPLACE, LDF/RESAMP.

%         Author:  J. Jezek  14-Jun-2000
%         Copyright(c) 1999 by Polyx, Ltd.
%         $ Revision $  $ Date 12-Jul-2000 $
%                       $ Date 16-Oct-2000 $
%                       $ Date 16-Jan-2001 $
%                       $ Date 03-Jan-2002 $
%                       $ Date 30-Sep-2002 $
%                       $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin<1,
   error('Not enough input arguments.');
end;
if ~isa(F,'ldf'),
   error('Some argument but not 1st is invalidly fraction.');
end;

T = 1; tau = 0; tol = PGLOBAL.ZEROING; readtau = 0;
met = 'ss'; var = PGLOBAL.DISCRVAR;

lv = length(varargin);
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
         if strcmp(arg,'zi'), arg = 'z^-1';
         end;
         I = strmatch(arg,{'z';'q';'z^-1';'d'},'exact');
         J = strmatch(arg,{'ss'},'exact');
         if ~isempty(I),
            var = arg;
         elseif ~isempty(J), 
            met = arg;
         else
            error(['Invalid ',nth(i+1),' argument.']);      
         end;
      end;
   elseif isa(arg,'double'),
      if ~isempty(arg),
         if ~readtau,
            if length(arg)==1 & isreal(arg) & arg>=0,
               T = arg; readtau = 1;
            else
               error(['Invalid ',nth(i+1),' argument.']);
            end;
         else
            if ndims(arg)==2 & any(size(arg)==1) & ...
                  isreal(arg) & all(arg>=0),
               tau = arg;
            else
               error(['Invalid ',nth(i+1),' argument.']);
            end;
         end;
      end;
      readtau = 1;
   else
      error(['Invalid ',nth(i+1),' argument.']);
   end;
end;

ltau = length(tau);
if isempty(F) | all(all(F.frac.num==0)),
   if ltau==1,
      G = F;
   else
      G = cell(size(tau));
      for i = 1:ltau,
         G{i} = F;
      end;
   end;
   return;
end;

Fvar = F.frac.v;
if ~isempty(Fvar) & ~strcmp(Fvar,'s') & ~strcmp(Fvar,'p'),
   error('Invalid variable symbol; must be continuous-time.');
end;

if strcmp(met,'ss'),
   [A,B,C,D] = abcd(F,tol);
   if D~=0 & any(tau==0),
      error('In time domain, sampling of Dirac impulse.');
   end;
   AA = expm(A*T); CC = C;
   if ltau==1,
      BB = expm(A*tau)*B;
   else
      [sB1,sB2] = size(B); sB2min1 = sB2-1;
      BB = zeros(sB1,0);
      for i = 1:ltau,
         BB = [BB, expm(A*tau(i))*B];
      end;
   end;
   
   GG = ldf(AA,BB,CC,num2str(tol),'z');
   GG.frac.num = shift(GG.frac.num,1,'z');
   GG.frac.h = T;
   if T>0,
      GG.frac.num.h = T; GG.frac.den.h = T;
   else
      GG.frac.num.h = []; GG.frac.den.h = [];
   end;
   
   if strcmp(var,'z^-1'),
      GG = reverse(GG);
   elseif strcmp(var,'q'),
      GG.frac.v = 'q';
   elseif strcmp(var,'d'),
      GG = reverse(GG); GG.frac.v = 'd';
   end;
   
   if strcmp(PGLOBAL.COPRIME,'cop'), GG = coprime(GG,tol);
   end;
   if strcmp(PGLOBAL.REDUCE,'red'), GG = reduce(GG,tol);
   else GG = smreduce(GG);
   end;
   if strcmp(PGLOBAL.DEFRACT,'defr'), GG = defract(GG);
   end;
   
   if ltau==1,
      G = GG;
   else
      G = cell(size(tau));
      if isa(GG,'ldf'),
         ii = 1; Gd = GG.frac.den;
         for i = 1:ltau,
            Gin = GG.frac.num(:,ii:ii+sB2min1);
            G{i} = ldf(Gd,Gin);
            ii = ii+sB2;
         end;
      else
         ii = 1;
         for i = 1:ltau,
            G{i} = GG(:,ii:ii+sB2min1);
            ii = ii+sB2;
         end;
      end;
   end;
   
end;
   
%end .. @ldf/samp
   