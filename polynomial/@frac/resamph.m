function G = resamph(F,varargin)
%RESAMPH     Resample and hold discrete-time fraction
%
% For fraction F with variable symbol 'z','q','z^-1' or 'd',
% the command  G = RESAMPH(F,K,H,L) returns fraction G with
% the same variable symbol, the result of holding with
% interval L and resampling with period K and phase H. Note
% that the meaning of the variable has changed:
%      new var = old var to the K-th power
%
% "Holding" or "zero order holding" with interval L means:
% In frequency domain: F is multiplied by
%     U(z) = (1-z^-L)/(1-z^-1)     
% In time-domain: F is convolved by U, which is
% held to 1 during the interval:
%     u(t) = 1  for t = 0,1,...L-1
%            0  for t = L,L+1,...
%
% The interpretation: F is the transfer function of a dynamic 
% system (in the old variable), G is the transfer function of
% the series:
%  -> holding element(L) -> system -> resampling element(K,H) ->
% (in the new variable).
%
% H can also be, instead of scalar, a vector. In such
% a case, the result is a cell vector.
%
% Arguments K,H,L are optional, the defaults are K = 2,
%  H = 0, L = K.
%
% For more details and for possible further arguments,
% see RDF/RESAMP, LDF/RESAMP, MDF/RESAMP, SDF/RESAMP.

%    Author: J.Jezek, 27-Jun-2000
%    Copyright(c) 2000 by Polyx,Ltd.
%    $ Revision $  $ Date 22-Aug-2000 $
%                  $ Date 04-Oct-2000 $
%                  $ Date 01-Dec-2000 $
%                  $ Date 01-Feb-2001 $
%                  $ Date 25-Jul-2002 $
%                  $ Date 30-Sep-2002 $
%                  $ Date 14-Oct-2002 $

global PGLOBAL;

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
   F = F*zz/(zz-1);
case 'q', qq = q(Fh);
   F = F*qq/(qq-1);
case 'z^-1', zizi = zi(Fh);
   F = F/(1-zizi);
case 'd', dd = d(Fh);
   F = F/(1-dd);
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
      case 'z', zz = z(Gh); 
         GG = GG*(zz-1)/zz;
      case 'q', qq = q(Gh);
         GG = GG*(qq-1)/qq;
      case 'z^-1', zizi = zi(Gh);
         GG = GG*(1-zizi);
      case 'd', dd = d(Gh);
         GG = GG*(1-dd);
      end;
      
      if lh==1, G = GG;
      else G{x} = GG;
      end;
   end;
   
else
   h = reshape(h,1,lh);
   hr2 = zeros(1,lh); quot2 = hr2;
   lk = floor(l/k); r2 = l-lk*k;
   for x = 1:lh,
      hr2x = h(x)-r2; quot2x = lk;
      if hr2x<0, hr2x = hr2x+k; quot2x = quot2x+1;
      end;
      hr2(x) = hr2x; quot2(x) = quot2x;
   end;
   
   hr = [h,hr2];
   eval('Ghr = resamp(F,k,hr,varargin{start:lv});', ...
      'error(peel(lasterr));');
   Ghr = reshape(Ghr,lh,2);
   Gv = Ghr{1}.v;
   if isempty(Gv), Gv = PGLOBAL.DISCRVAR;
   end;
   Gh = Ghr{1}.h;
   
   for x = 1:lh,
      quot2x = quot2(x);
      switch Gv,
      case 'z^-1', zizi = zi(Gh);
         GG = Ghr{x,1} - zizi^quot2x*Ghr{x,2};
      case 'd', dd = d(Gh);
         GG = Ghr{x,1} - dd^quot2x*Ghr{x,2};
      case 'z', zz = z(Gh); zzquot = zz^quot2x;
         GG = (zzquot*Ghr{x,1} - Ghr{x,2})/zzquot;
      case 'q', qq = q(Gh); qqquot = qq^quot2x;
         GG = (qqquot*Ghr{x,1} - Ghr{x,2})/qqquot;      
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

%end .. @frac/resamph

