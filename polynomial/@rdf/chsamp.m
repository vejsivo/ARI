function H = chsamp(G,varargin)
%CHSAMP     Change sampling of discrete time
%                  right-den fraction
%
% Discrete time right-den fraction G is supposed to
% be a result of sampling of continuous time fraction
% F with period T and phase TAU. The command
%    H = CHSAMP(G,T,TN,TAU,TAUN)   returns the result
% of sampling F with new period TN and new phase TAUN.
%
% Fraction G, as a function of complex variable,
% must be proper, i.e. it must not have a pole for
% z=Inf (z^-1=0, q=Inf, d=0). Furthermore, G must have
% a zero point for z=0 (z^-1=Inf, q=0, d=Inf). The only
% case when these conditions are not required, is
% TN/T integer & TAU=0 & TAUN/T integer.
%
% The result H is function of the same variable as G,
% formally. However, note that the meaning of the
% variable has changed:
%     new var = old var to the TN/T -th power
%
% When TN/T is not integer, the changing of sampling
% is not unique. For the result H, the main value is
% selected, i.e. the complex function whose poles
% in 'z' have argument (the angle in polar coordinates)
% closest to zero.
% 
% For real fraction G, the resulting fraction H is
% usually also real. However, when G has a real
% negative pole in 'z' and TN/T is not integer, the
% result may happen to be complex.
%
% The arguments T,TN,TAU,TAUN are optional, they may
% be omitted of gien by []. Their defaults are
%  T = G.h, TN = T/2, TAU = 0, TAUN = 0.
% The argument TAUN can also be, instead of scalar,
% a vector. In such a case, the result is a cell vector
% whose entries contain results for individual TAUN's.
% All cells have the same denominator.
%
% In case of zero sampling period, only special G is
% allowed: it must be equal to  K*z/(z-1)  or
%  K*q/(q-1),  K/(1-z^-1),  K/(1-d) . The result H
% is equal to G.
%
% An optional argument TOL may specify
% the zeroing tolerance to be used instead of the
% standard one. It must have a form of character
% string, e.g. '10^-6'. An optional argument MET
% may specify the method used. The methods are
%    'ss' - state space
%
% See also RDF/SAMP, FRAC/UNSAMP, RDF/RESAMP.

%      Author: J.Jezek, 05-Dec-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 05-Feb-2001 $
%                    $ Date 14-Oct-2002 $
%                    $ Date 28-Feb-2003 $

global PGLOBAL;

if nargin<1,
   error('Not enough input arguments.');
end;
if ~isa(G,'rdf'),
   error('Some argument but not 1st is invalidly fraction.');
end;

lv = length(varargin); 
T = G.frac.h; Tn = T/2; tau = 0; taun = 0;
start = 1;
if lv>=1,
   arg = varargin{1};
   if isa(arg,'double') & ndims(arg)==2,
      if ~isempty(arg),
         if length(arg)==1 & isreal(arg) & arg>=0,
            T = arg; Tn = T/2;
         else
            error('Invalid old sampling period.');
         end;
      end;
      start = 2;
      if lv>=2,
         arg = varargin{2};
         if isa(arg,'double') & ndims(arg)==2,
            if ~isempty(arg),
               if length(arg)==1 & isreal(arg) & arg>=0,
                  Tn = arg;
               else
                  error('Invalid new sampling period.');
               end;
            end;
         end;
         start = 3;
         if lv>=3,
            arg = varargin{3};
            if isa(arg,'double') & ndims(arg)==2,
               if ~isempty(arg),
                  if length(arg)==1 & isreal(arg) & arg>=0,
                     tau = arg;
                  else
                     error('Invalid old sampling phase.');
                  end;
               end;
               start = 4;
               if lv>=4,
                  arg = varargin{4};
                  if isa(arg,'double') & ndims(arg)==2,
                     if ~isempty(arg),
                        if any(size(arg)==1) & ...
                              isreal(arg) & all(arg>=0),
                           taun = arg;
                        else
                           error('Invalid new sampling phase.');
                        end;
                     end;
                     start = 5;
                  end;
               end;
            end;
         end;
      end;
   end;
end;

met = 'ss'; tol = PGLOBAL.ZEROING;
for i = start:lv,
   arg = varargin{i};
   if isa(arg,'char'),
      num = str2num(arg);
      if isa(num,'double') & ~isempty(num),
         if length(num)==1 & isreal(num) & num>=0 & num<=1,
            tol = num;
         else
            error(['Invalid ',nth(i+1),' argument.']);
         end;
      elseif ~isempty(arg),
         met = arg;
      end;
   else
      error(['Invalid ',nth(i+1),' argument.']);
   end;
end;
strtol = num2str(tol);

ltaun = length(taun); staun = size(taun);
HH = G;
if ltaun==1,
   H = HH;
else
   H = cell(staun);
   for k = 1:ltaun,
      H{k} = HH;
   end;
end;
if isempty(G) | all(all(G.frac.num==0)),
   return;
end;

Gvar = G.frac.v;
if strcmp(Gvar,'s') | strcmp(Gvar,'p'),
   error('Invalid variable symbol in 1st argument; must be discrete time.');
end;

if isempty(T),
   error('Old sampling period missing.');
end;
if ~isfinite(T),
   error('Old sampling period is not finite.');
end;

if T>0,
   kT = Tn/T; ktau = taun/T;
   G.frac.h = T;
   if kT>0 & floor(kT)==kT & tau==0 & all(floor(ktau)==ktau),
      eval('H = resamp(G,kT,ktau,varargin{start:lv});', ...
         'error(peel(lasterr));');
      return;
   end;
end;

if ~isproper(G),
   error('Invalid 1st argument; is not proper.');
end;   

switch Gvar,
case 'z', text = 'z=0';
case 'q', text = 'q=0'; G.frac.v = 'z';
case 'z^-1', text = 'z^-1=inf'; G = reverse(G);
case 'd', text = 'd=inf'; G.frac.v = 'z^-1'; G = reverse(G);
case '',
   error('Invalid 1st argument; has not required zero point.');
end;
absterm = G.frac.num{0};
if any(any(absterm~=0)),
   error(['Invalid 1st argument; has not zero point ',text,'.']);
end;

if T==0,
   Gh = G.frac.h;
   if isempty(Gh) | Gh~=0,
      num = shift(G.frac.num,-1);
      if isa(num,'tsp'),
         error('Invalid 1st argument for zero old sampling period.');
      end;
      zz = z(num.h);
      num = num*(zz-1); den = G.frac.den;
      [Q,R] = rdiv(num,den);
      if deg(R)>=0 | deg(Q)>0,
         error('Invalid 1st argument for zero old sampling period.');
      end;
   end;
   G.frac.h = Tn;
   if ltaun==1,
      H = G;
   else
      for k = 1:ltaun,
         H{k} = G;
      end;
   end;
   return;
end;

if Tn==0,
   eval(['F = unsamp(G,T,tau,varargin{start:lv});', ...
         'H = samp(F,Tn,taun,varargin{start:lv});'], ...
      'error(peel(lasterr));');
   return;
end;

if strcmp(met,'ss'),
   G.frac.num = shift(G.frac.num,-1);
   [AA,BB,CC,DD] = abcd(G,tol);
   AAn = AA^kT;
   BBn = BB;
   if ltaun==1,
      if taun==tau, CCn = CC;
      else CCn = CC*AA^((taun-tau)/T);
      end;
   else
      [sCC1,sCC2] = size(CC); sCC1min1 = sCC1-1;
      CCn = zeros(0,sCC2);
      for i = 1:ltaun,
         CCn = [CCn; CC*AA^((taun(i)-tau)/T)];
      end;
   end;
   
   HH = rdf(AAn,BBn,CCn,strtol,'z');
   HH.frac.num = shift(HH.frac.num,1,'z');
   HH.frac.h = Tn;
   HH.frac.num.h = Tn; HH.frac.den.h = Tn;
   
   if strcmp(Gvar,'z^-1'),
      HH = reverse(HH);
   elseif strcmp(Gvar,'q'),
      HH.frac.v = 'q';
   elseif strcmp(Gvar,'d'),
      HH = reverse(HH); HH.frac.v = 'd';
   end;
   
   if strcmp(PGLOBAL.COPRIME,'cop'), HH = coprime(HH,tol);
   end;
   if strcmp(PGLOBAL.REDUCE,'red'), HH = reduce(HH,tol);
   else HH = smreduce(HH);
   end;
   if strcmp(PGLOBAL.DEFRACT,'defr'), HH = defract(HH);
   end;
   
   if ltaun==1,
      H = HH;
   else
      H = cell(staun);
      ii = 1; Hd = HH.frac.den;
      for i = 1:ltaun,
         Hin = HH.frac.num(ii:ii+sCC1min1,:);
         H{i} = rdf(Hin,Hd);
         ii = ii+sCC1;
      end;
   end;
   
else
   error('Invalid command option.');
end;
 
%end .. @rdf/chsamp
 
