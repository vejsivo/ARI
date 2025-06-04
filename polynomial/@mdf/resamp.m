function G = resamp(F,varargin)
%RESAMP    Resample discrete time matrix-den fraction
%
% Let F be a matrix-denominator fraction in variable 'z','q','z^-1'
% or 'd'; its Laurent series in the point  z=Inf  being
%   F(-m)*z^m + ... + F(0) + ... + F(n)*z^-n + ...
% The command  G = RESAMP(F,K,H)  where scalar integers K,H are
% resampling period, K >= 1, and resampling phase, H >= 0,
% (defaults being  K=2, H=0), returns matrix-denominator fraction G
% in the same variable whose Laurent series equals to the resampled
% series of F. Note that the meaning of the variable has changed:
%    new var = old var  to the K-th power
% The sampling period  G.h  of G is K-times greater than  F.h .
%
% For example:
%    F = z./(z-0.8) = 1./(1-0.8*z^-1)
%    laurent = 1 + 0.8*z^-1 + 0.64*z^-2 + 
%              0.512*z^-3 + 0.0496*z^-4 + 0.32768*z^-5 + ...
%
%    RESAMP(F,3,0) = z./(z-0.512) = 1./(1-0.512*z^-1)
%    laurent = 1 + 0.512*z^-1 + ...
%
%    RESAMP(F,3,1) = (0.8*z)./(z-0.512) = 0.8./(1-0.512*z^-1)
%    laurent = 0.8 + 0.4096*z^-1 + ...
%
%    RESAMP(F,3,2) = (0.64*z)./(z-0.512) = 0.64./(1-0.512*z^-1)
%    laurent = 0.64 + 0.32768*z^-1 + ...
%
% H can also be, instead of scalar, a vector. In this case,
% the results form a cell vector.
%
% An optional argument TOL may specify a zeroing tolerance to be
% used instead of the standard one. It must have a form of character
% string, e.g. '10^-6' .
%
% An optional argument MET may specify the method to be used.
% The methods are 'ss' - state space,
%                 'cr1' - complex roots of one.
%
% See also MDF/LAURENT, POL/RESAMP, TSP/RESAMP.

%         Author:  J. Jezek  10-Dec-1999
%         Copyright(c) 1999 by Polyx, Ltd.
%         $ Revision $  $ Date 25-Apr-2000 $
%                       $ Date 25-Jul-2000 $
%                       $ Date 16-Oct-2000 $
%                       $ Date 01-Feb-2001 $
%                       $ Date 14-Oct-2002 $
%                       $ Date 28-Feb-2003 $

global PGLOBAL;

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ~isa(F,'mdf'),
   error('Some argument but not 1st is invalidly fraction.');
end;

k = 2; h = 0; tol = PGLOBAL.ZEROING; readh = 0;
met = 'ss';

lv = length(varargin);
for i = 1:lv,
   arg = varargin{i};
   if isa(arg,'char'),
      num = str2num(arg);
      if isa(num,'double') & ~isempty(num),
         if length(num)==1 & isreal(num) & num>=0 & num<=1,
            tol = num;
         else
            error('Invalid tolerance.');
         end;
      elseif ~isempty(arg),
         met = arg;
      end;
   elseif isa(arg,'double'),
      if ~isempty(arg),
         if ~readh,
            if length(arg)==1 & isreal(arg) & arg>0,
               k = arg;
            else
               error('Invalid resampling period.');
            end;
         else
            if ndims(arg)==2 & any(size(arg)==1) & ...
                  isreal(arg) & all(arg>=0),
               h = arg;
            else
               error('Invalid resampling phase.');
            end;
         end;
      end;
      readh = 1;
   else
      error(['Invalid ',nth(i+1),' argument.']);
   end;
end;

if k~=floor(k),
   error('Invalid resampling period; must be integer.');
end;
if h~=floor(h),
   error('Invalid resampling phase; must be integer.');
end;

var = F.frac.v;
if ~isempty(var),
   I = strmatch(var,{'z';'z^-1';'q';'d'},'exact');
   if isempty(I),
      error('Invalid variable symbol; must be discrete time.');
   end;
end;

if ~isempty(F.frac.h), Gh = F.frac.h*k;
else Gh = [];
end;

pvar = pol([0,1],1,var); pvar.h = F.frac.h;
switch var,
case {'z','q'}, avar = 'z^-1';
otherwise avar = 'z';
end;

[sF1,sF2] = size(F);
Q = pol(zeros(sF1,sF2)); X = Q;
t = tdeg(F.frac.den,'ent');
for u = 1:sF1,
   for v = 1:sF2,
      Fd = F.frac.den(u,v);
      [Q(u,v),Fn] = rdiv(F.frac.num(u,v),Fd,tol);
      tt = t(u,v);
      Fd = shift(Fd,-tt);
      F.frac.den(u,v) = Fd;
      [XX,F.frac.num(u,v)] = axbyc(Fd,pvar^tt,Fn,'miny');
      X(u,v) = rev(XX,tt);
   end;
end;
X.v = avar;

lh = length(h);
if strcmp(met,'cr1'),
   
   kf = factor(k); lkf = length(kf);
   sh = size(h); G = cell(sh);
   for y = 1:lh, G{y} = F;
   end;
   hx = h;
   for x = 1:lkf,
      kfx = kf(x); hfx = mod(hx,kfx);
      
      if kfx==2,
         for y = 1:lh,
            Gy = G{y};
            bb = Gy.frac.num; aa = Gy.frac.den;
            dbb = deg(bb); daa = deg(aa);
            for i = 1:2:dbb,
               bb{i} = -bb{i};
            end;
            for i = 1:2:daa,
               aa{i} = -aa{i};
            end;
            Gyalt = mdf(bb,aa);
            if hfx(y)==0, H = plus(Gy,Gyalt,tol);
            else H = minus(Gy,Gyalt,tol);
            end;
            G{y} = H*0.5;
         end;
         
      else 
         omega = exp(2*pi*j/kfx);
         if strcmp(var,'z^-1') | strcmp(var,'d'),
            omegatran = omega;
         else
            omegatran = 1/omega;
         end;
   
         k1 = kfx-1; H = cell(sh);
         for y = 1:lh,
            omegakh = omega^(kfx-hfx(y));
            H{y} = zeros(F.frac.s);
            for i = 0:k1,
               omegai = omegatran^i;
               bb = linvt(G{y}.frac.num,omegai,0);
               aa = linvt(G{y}.frac.den,omegai,0);
               H{y} = plus(H{y}, omegakh^i * mdf(bb,aa), tol);
            end;
            G{y} = H{y}/kfx;
         end;
      end;
      
      if strcmp(var,'z^-1') | strcmp(var,'d'),
         for y = 1:lh,
            Gy = G{y}; Gyn = Gy.frac.num; Gyd = Gy.frac.den;
            oo = min(tdeg(Gyn),tdeg(Gyd)+hfx(y));
            if ~isempty(oo),
               Gyn = shift(Gyn,-oo); Gyd = shift(Gyd,hfx(y)-oo);
            end;
            G{y} = mdf(resamp(Gyn,kfx,0),resamp(Gyd,kfx,0));
         end;
      else
         for y = 1:lh,
            Gy = G{y}; Gyn = Gy.frac.num; Gyd = Gy.frac.den;
            oo = min(tdeg(Gyn)+hfx(y),tdeg(Gyd));
            if ~isempty(oo),
               Gyn = shift(Gyn,hfx(y)-oo); Gyd = shift(Gyd,-oo);
            end;
            G{y} = mdf(resamp(Gyn,kfx,0),resamp(Gyd,kfx,0));
         end;
      end;
         
      hx = (hx-hfx)/kfx;
   end;   %for x
   
   for y = 1:lh,
      G{y} = G{y} + resamp(Q,k,h(y));
   end;
         
   Xd = deg(X);
   if ~isempty(Xd) & Xd>=0,
      pvarXd = pvar^Xd;
      for y = 1:lh,
         Xs = resamp(X,k,h(y));
         Xs = rev(Xs,Xd);
         pvarXd.h = Xs.h;
         G{y} = G{y} + mdf(Xs,pvarXd);
      end;
   end;
   
   if lh==1, G = G{1};
   end;
      
elseif strcmp(met,'ss'),

   [A,B,C,D] = abcd(F,tol);
   AA = A^k; BB = A^(k-1)*B;
   if lh == 1,
      CC = C*A^h;
      if h==0, DD = D;
      else DD = C*A^(h-1)*B;
      end;
   else
      [sC1,sC2] = size(C); [sD1,sD2] = size(D);
      CC = zeros(0,sC2); DD = zeros(0,sD2);
      for i = 1:lh;
         hi = h(i);
         if hi==0, CC = [CC;C]; DD = [DD;D];
         else Ahi1 = A^(hi-1);
            CC = [CC; C*A*Ahi1];
            DD = [DD; C*Ahi1*B];
         end;
      end;
   end;
   
   GG = mdf(AA,BB,CC,DD,num2str(tol),'z');
   GG.frac.h = Gh;
   
   if strcmp(var,'z^-1'),
      GG = reverse(GG);
   elseif strcmp(var,'q'),
      GG.frac.v = 'q';
   elseif strcmp(var,'d'),
      GG = reverse(GG); GG.frac.v = 'd';
   end;
   
   [sQ1,sQ2] = size(Q); sQ1min1 = sQ1-1;
   QQ = pol(zeros(0,sQ2));
   for y = 1:lh,
      QQ = [QQ; resamp(Q,k,h(y))];
   end;
   GG = GG + QQ;

   Xd = deg(X);
   if ~isempty(Xd) & Xd>=0,
      XX = pol(zeros(0,sQ2));
      for y = 1:lh,
         Xs = resamp(X,k,h(y));
         Xs = rev(Xs,Xd); Xs.v = var;
         XX = [XX;Xs];
      end;
      pvar.h = XX.h;
      GG = GG + mdf(XX,pvar^Xd);
   end;

   if lh==1,
      G = GG;
   else
      G = cell(size(h));
      ii = 1;
      for i = 1:lh,
         if isa(GG,'mdf'),
            Gin = GG.frac.num(ii:ii+sQ1min1,:);
            Gid = GG.frac.den(ii:ii+sQ1min1,:);
            G{i} = Gin./Gid;
         else  % isa(GG,'pol')
            G{i} = GG(ii:ii+sQ1min1,:);
         end;
         ii = ii+sQ1;
      end;
   end;
   
else
   error('Invalid command option.');
end;
   
%end .. @mdf/resamp
