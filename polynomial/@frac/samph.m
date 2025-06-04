function G = samph(F,varargin)
%SAMPH    Sample and hold continuous time fraction
%
% For fraction F in variable 's' or 'p', the command
%  G = SAMPH(F,T,TAU,L) with T,L > 0, TAU >= 0, returns
% proper fraction G in the variable 'z','q','z^-1' or 'd',
% the result of holding with interval L and sampling with
% period T and phase TAU.
%
% "Holding" or "zero order holding" with interval L means:
% In frequency domain: F is multiplied by
%     U(s) = (1-exp(-s*L))/s
% In time domain: F is convolved by U, which is
% held to 1 during the interval (0,L) :
%     u(t) = 1  for 0 <= t < L
%            0  for L <= t
%
% The interpretation: F is the transfer function of a dynamic 
% system (continuous time), G is the transfer function of the
% series:
%  -> holding element(L) -> system -> sampling element(T,TAU) ->
% (discrete time).
%
% TAU can also be, instead of scalar, a vector. In this case,
% the results form a cell vector. All cells have the
% same denominator.
%
% Arguments T,TAU,L are optional, the defaults being T=1,TAU=0,
% L=T. When some of the sampling points TAU, T+TAU, 2*T+TAU ...
% coincides with 0 or L, the fraction F must be proper (otherwise,
% in time domain, sampling of a Dirac impulse would be attempted).
%
% For more details and for possible further arguments,
% see RDF/SAMP, LDF/SAMP, MDF/SAMP, SDF/SAMP.

%    Author: J.Jezek, 13-Jul-2000
%    Copyright(c) 2000 by Polyx,Ltd.
%    $ Revision $  $ Date 04-Oct-2000 $
%                  $ Date 01-Dec-2000 $
%                  $ Date 20-Jan-2001 $
%                  $ Date 30-Sep-2002 $
%                  $ Date 14-Oct-2002 $
%                  $ Date 28-Feb-2003 $

global PGLOBAL;

if nargin<1,
   error('Not enough input arguments.');
end;
if ~isa(F,'frac'),
   error('Some argument but not 1st is invalidly fraction.');
end;

T = 1; tau = 0; L = T;
lv = length(varargin); start = 1;
if lv>=1,
   arg1 = varargin{1};
   if isa(arg1,'double'),
      if ~isempty(arg1),
         if length(arg1)==1 & isreal(arg1) & arg1>0,
            T = arg1; L = T;
         else
            error('Invalid sampling period.');
         end;
      end;
      start = 2;
      if lv>=2,
         arg2 = varargin{2};
         if isa(arg2,'double'),
            if ~isempty(arg2),
               if ndims(arg2)==2 & any(size(arg2)==1) & ...
                     isreal(arg2) & all(arg2>=0),
                  tau = arg2;
               else
                  error('Invalid sampling phase.');
               end;
            end;
            start = 3;
            if lv>=3,
               arg3 = varargin{3};
               if isa(arg3,'double'),
                  if ~isempty(arg3),
                     if length(arg3)==1 & isreal(arg3) & arg3>0,
                        L = arg3;
                     else
                        error('Invalid holding interval.');
                     end;
                  end;
                  start = 4;
               end;
            end;
         end;
      end;
   end;
end;

tol = PGLOBAL.ZEROING;
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
      end;
   end;
end;

ltau = length(tau); stau = size(tau);
GG = F;
if ltau==1,
   G = GG;
else
   G = cell(stau);
   for k = 1:ltau,
      G{k} = GG;
   end;
end;
if isempty(F) | all(all(F.num==0)),
   return;
end;

sdefr = PGLOBAL.DEFRACT; PGLOBAL.DEFRACT = 'ndefr';
scop = PGLOBAL.COPRIME; PGLOBAL.COPRIME = 'ncop';
Fv = F.v;
if isempty(Fv), Fv = PGLOBAL.CONTVAR;
end;

switch Fv,
case 's',
   F.den = F.den*s;
case 'p',
   F.den = F.den*p;
otherwise,
   error('Invalid variable symbol; must be continuous-time.');
end;

if L==T,
   eval('G = samp(F,T,tau,varargin{start:lv});', ...
         ['PGLOBAL.DEFRACT = sdefr;', ...
         'PGLOBAL.COPRIME = scop;', ...
         'error(peel(lasterr));']);
   for k = 1:ltau,
      if ltau==1, GG = G;
      else GG = G{k};
      end;
      Gv = GG.v;
      if isempty(Gv), Gv = PGLOBAL.DISCRVAR;
      end;
      Gh = GG.h; Gn = GG.num; Gd = GG.den;
         
      switch Gv,
      case 'z^-1', zizi = zi(Gh);
         Gn = Gn*(1-zizi);
      case 'd', dd = d(Gh);
         Gn = Gn*(1-dd);
      case 'z', zz = z(Gh);
         Gn = Gn*(zz-1); Gd = Gd*zz;
      case 'q', qq = q(Gh);
         Gn = Gn*(qq-1); Gd = Gd*qq;
      end;
      
      GG.num = Gn; GG.den = Gd; GG.v = Gv; GG.h = Gh;
      if ltau==1, G = GG;
      else G{k} = GG;
      end;
   end;
   
else
   tau = reshape(tau,1,ltau);
   taur2 = zeros(1,ltau); quot2 = taur2;
   LT = floor(L/T); r2 = L-LT*T;
   for k = 1:ltau,
      taur2k = tau(k)-r2; quot2k = LT;
      if taur2k<0, taur2k = taur2k+T; quot2k = quot2k+1;
      end;
      taur2(k) = taur2k; quot2(k) = quot2k;
   end;
    
   taur = [tau, taur2];
   eval('Gtaur = samp(F,T,taur,varargin{start:lv});', ...
         ['PGLOBAL.DEFRACT = sdefr;', ...
         'PGLOBAL.COPRIME = scop;', ...
         'error(peel(lasterr));']);
   Gtaur = reshape(Gtaur,ltau,2);
   Gv = Gtaur{1}.v;
   if isempty(Gv), Gv = PGLOBAL.DISCRVAR;
   end;
   Gh = Gtaur{1}.h;
   
   for k = 1:ltau,
      quot2k = quot2(k);
      switch Gv,
      case 'z^-1', zizi = zi(Gh);
         Gn = Gtaur{k,1}.num - zizi^quot2k*Gtaur{k,2}.num;
         Gd = Gtaur{k,1}.den;
      case 'd', dd = d(Gh);
         Gn = Gtaur{k,1}.num - dd^quot2k*Gtaur{k,2}.num;
         Gd = Gtaur{k,1}.den;
      case 'z', zz = z(Gh); zzquot = zz^quot2k;
         Gn = zzquot*Gtaur{k,1}.num - Gtaur{k,2}.num;
         Gd = zzquot*Gtaur{k,1}.den;
      case 'q', qq = q(Gh); qqquot = qq^quot2k;
         Gn = qqquot*Gtaur{k,1}.num - Gtaur{k,2}.num;
         Gd = qqquot*Gtaur{k,1}.den;
      end;
      
      GG.num = Gn; GG.den = Gd; GG.v = Gv; GG.h = Gh;
      if ltau==1, G = GG;
      else G{k} = GG;
      end;
   end;         
end;

PGLOBAL.COPRIME = scop; PGLOBAL.DEFRACT = sdefr;
for k = 1:ltau,
  if ~isa(G,'cell'), GG = G;      
  else GG = G{k};
  end;
  GG.c = 'cop?'; GG.r = 'red?';
  if strcmp(PGLOBAL.COPRIME,'cop'),
     GG = coprime(GG,tol);
  end;
  if strcmp(PGLOBAL.REDUCE,'red'),
     GG = reduce(GG,tol);
  else
     GG = smreduce(GG);
  end;
  if strcmp(PGLOBAL.DEFRACT,'defr'),
     GG = defract(GG);
  end;
  if ~isa(G,'cell'), G = GG;
  else G{k} = GG;
  end;
end;
      
%end .. @frac/samph
