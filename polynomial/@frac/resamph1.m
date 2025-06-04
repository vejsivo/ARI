function G = resamph1(F,varargin)
%RESAMPH1    Resample and first order hold discrete time fraction
%
% For fraction F in variable 'z','q','z^-1' or 'd', the
% command  G = RESAMPH1(F,K,H,L) returns fraction G in the
% same variable, the result of first order holding with
% interval L and resampling with period K and phase H. Note
% that the meaning of the variable has changed:
%    new var = old var to the K-th power
%
% "First order holding" with interval L means:
% In frequency domain: F is multiplied by
%     U(z) = z^-1*(1-z^-L)^2/(1-z^-1)^2
% In time-domain: F is convolved by U, which is the first
% order spline:
%  u(t) =  t                            for t = 0,1,...L-1
%          t - 2*(t-L)                  for t = L,L+1,...2*L-1
%          t - 2*(t-L) + (t-2*L) = 0    for t = 2*L,2*L+1,...
%
% The interpretation: F is the transfer function of a dynamic 
% system (in the old variable), G is the transfer function of
% the series:
%  -> first order holding element(L) -> system ->
%         -> resampling element(K,H) ->
% (in the new variable).
%
% Arguments K,H,L are optional, the defaults are K = 2,
%  H = 0, L = K.
%
% For more details and for possible further arguments,
% see RDF/RESAMP, LDF/RESAMP, MDF/RESAMP, SDF/RESAMP,
% FRAC/RESAMPH.

%    Author: J.Jezek, 27-Jun-2000
%    Copyright(c) 2000 by Polyx,Ltd.
%    $ Revision $  $ Date 22-Aug-2000 $
%                  $ Date 04-Oct-2000 $
%                  $ Date 01-Dec-2000 $
%                  $ Date 02-Feb-2001 $
%                  $ Date 25-Jul-2002 $
%                  $ Date 30-Sep-2002 $
%                  $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin<1,
   error('Not enough input arguments.');
end;
Fcl = class(F);
if strcmp(Fcl,'frac'),
   error('Invalid 1st argument.');
end;
if ~isa(F,'frac'),
   error('Some argument but not 1st is invalidly fraction.');
end;

k = 2; h = 0; l = k;
lv = length(varargin); start = 1;
if lv>=1,
   arg1 = varargin{1};
   if isa(arg1,'double'),
      if ~isempty(arg1),
         if length(arg1)==1 & isreal(arg1) & ...
               arg1>0 & floor(arg1)==arg1,
            k = arg1; l = k;
         else
            error('Invalid resampling period.');
         end;
      end;
      start = 2;
      if lv>=2,
         arg2 = varargin{2};
         if isa(arg2,'double'),
            if ~isempty(arg2),
               if ndims(arg2)==2 & any(size(arg2)==1) & ...
                     isreal(arg2) & all(floor(arg2)==arg2) & ...
                     all(arg2>=0),
                  h = arg2;
               else
                  error('Invalid resampling phase.');
               end;
            end;
            start = 3;
            if lv>=3,
               arg3 = varargin{3};
               if isa(arg3,'double'),
                  if ~isempty(arg3),
                     if length(arg3)==1 & isreal(arg3) & ...
                           floor(arg3)==arg3 & arg3>=1,
                        l = arg3;
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

lh = length(h); sh = size(h);
GG = F;
if lh==1,
   G = GG;
else
   G = cell(sh);
   for x = 1:lh,
      G{x}= GG;
   end;
end;
if isempty(F) | all(all(F.num==0)),
   return;
end;

sdefr = PGLOBAL.DEFRACT; PGLOBAL.DEFRACT = 'ndefr';
scop = PGLOBAL.COPRIME; PGLOBAL.COPRIME = 'ncop';
Fv = F.v;
if isempty(Fv), Fv = PGLOBAL.DISCRVAR;
end;
Fh = F.h;
if isempty(Fh) | ~isfinite(Fh), Fh = 1;
end;

switch Fv,
case 'z', zz = z(Fh);
   F = F*zz/(zz-1)^2;
case 'q', qq = q(Fh);
   F = F*qq/(qq-1)^2;
case 'z^-1', zizi = zi(Fh);
   F = F*zizi/(1-zizi)^2;
case 'd', dd = d(Fh);
   F = F*dd/(1-dd)^2;
otherwise,
   error('Invalid variable symbol; must be discrete-time.');
end;
F.h = Fh;

if l==k,
   eval('G = resamp(F,k,h,varargin{start:lv});',...
      'error(peel(lasterr));');
   for x = 1:lh,
      if lh==1, GG = G;
      else GG = G{x};
      end;   
      Gv = GG.v;
      if isempty(Gv), Gv = PGLOBAL.DISCRVAR;
      end;
      Gh = GG.h;
   
      switch Gv,
      case 'z',    zz = z(Gh);      GG = GG*(zz-1)^2/zz^2;
      case 'q',    qq = q(Gh);      GG = GG*(qq-1)^2/qq^2;
      case 'z^-1', zizi = zi(Gh);   GG = GG*(1-zizi)^2;
      case 'd',    dd = d(Gh);      GG = GG*(1-dd)^2;
      end;
      
      if lh==1, G = GG;
      else G{x} = GG;
      end;
   end;
      
else
   h = reshape(h,1,lh);
   hr2 = zeros(1,lh); quot2 = hr2;
   hr3 = hr2; quot3 = quot2;
   lk = floor(l/k); r2 = l-lk*k;
   ll = 2*l; llk = floor(ll/k); r3 = ll-llk*k;
   for x = 1:lh,
      hx = h(x);
      hr2x = hx-r2; quot2x = lk;
      if hr2x<0, hr2x = hr2x+k; quot2x = quot2x+1;
      end;
      hr2(x) = hr2x; quot2(x) = quot2x;
      hr3x = hx-r3; quot3x = llk;
      if hr3x<0, hr3x = hr3x+k; quot3x = quot3x+1;
      end;
      hr3(x) = hr3x; quot3(x) = quot3x;
   end;
   
   hr = [h,hr2,hr3];
   eval('Ghr = resamp(F,k,hr,varargin{start:lv});', ...
      'error(peel(lasterr));');
   Ghr = reshape(Ghr,lh,3);
   Gv = Ghr{1}.v;
   if isempty(Gv), Gv = PGLOBAL.DISCRVAR;
   end;
   Gh = Ghr{1}.h;
   
   for x = 1:lh,
      quot2x = quot2(x); quot3x = quot3(x);
      switch Gv,
      case 'z^-1', zizi = zi(Gh);
         GG = Ghr{x,1} - 2*zizi^quot2x*Ghr{x,2} + zizi^quot3x*Ghr{x,3};
      case 'd', dd = d(Gh);
         GG = Ghr{x,1} - 2*dd^quot2x*Ghr{x,2} + dd^quot3x*Ghr{x,3};
      case 'z', zz = z(Gh);
         quot = max(quot2x,quot3x); zzquot = zz^quot;
         quot2x = quot-quot2x; quot3x = quot-quot3x;
         GG = (zzquot*Ghr{x,1} - 2*zz^quot2x*Ghr{x,2} + zz^quot3x*Ghr{x,3})/zzquot;
      case 'q', qq = q(Gh);
         quot = max(quot2x,quot3x); qqquot = qq^quot;
         quot2x = quot-quot2x; quot3x = quot-quot3x;
         GG = (qqquot*Ghr{x,1} - 2*qq^quot2x*Ghr{x,2} + qq^quot3x*Ghr{x,3})/qqquot;
      end;
      
      GG.v = Gv; GG.h = Gh;
      if lh==1, G = GG;
      else G{x} = GG;
      end;
   end;  
end;

PGLOBAL.DEFRACT = sdefr; PGLOBAL.COPRIME = scop;
for x = 1:lh,
   if ~isa(G,'cell'), GG = G;
   else GG = G{x};
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
   else G{x} = GG;
   end;
end;

%end .. @frac/resamph1
