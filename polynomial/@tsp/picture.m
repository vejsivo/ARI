function H = picture(A,degree,T,tau,varargin)
%PICTURE    Picture two-sided polynomial
%
% For scalar two-sided polynomial A, the command
%  PICTURE(A)  pictures A, i.e. plots its
% coefficients  A{K}  versus K, in powers of z^-1
% (to represent time sigals in two-sided Z-transform).
% The command  PICTURE(A,DEGREE,T,TAU)  plots the
% coefficients in range  -DEGREE:DEGREE  versus
% time  K*T + TAU . If DEGREE is vector
%  [DEG1,DEG2]  then the range of K is  DEG1:DEG2 .
% Arguments DEGREE,T,TAU may be omitted or given
% by []. Default for DEGREE is  [TDEG(A),DEG(A)]
% where A is considered in powers of z^-1. Default
% for T is  A.h  or 1, for TAU is 0. Argument TAU
% can also be, instead of scalar, a vector; in such
% a case, for all TAUs, the same value is pictured.
%
% For matrix two-sided polynomial A, the command
% pictures all entries  A(I,J)  as said above,
% each entry in a separate figure.
%
% Standardly, new figures are created and their
% handles H(i,j), natural numbers, may be returned
% in the output argument  H = PICTURE(A) . 
% However, when an optional input argument H,
% matrix of natural numbers, is used in
%  PICTURE(A,T,TAU,DEGREE,H) , the existing
% figures may be used.
%
% By optional input argument S, character string,
% various line types, plot symbols and colors may
% be specified, for all entries of A uniformly,
% as in PLOT. The default is 'bs' (blue square)
% when sampling period is nonzero, and 'ks'
% (black square) in other cases.
%
% If A is complex then the real and the imaginary 
% part are pictured in the same figure, by default
% the former one by squares and the latter one by
% diamonds. Instead of optional input argument S,
% two arguments S1,S2 may be used.
%
% See also POL/PICTURE, PLOT, FIGURE.

%         Author:  J.Jezek, 02-Jan-2001
%         Copyright(c) 2001 by Polyx, Ltd.

ni = nargin; no = nargout;
if ni<1,
   error('Not enough input arguments.');
end;
if ~isa(A,'tsp'),
   error('Invalid 1st argument.');
end;

if isempty(A),
   if no, H = zeros(size(A));
   end;
   return;
end;

As = size(A); As1 = As(1); As2 = As(2);
Ah = A.h; Atau = 0;
Atdeg = -deg(A); Adeg = -tdeg(A);
Ad = Adeg; Atd = Atdeg;
HH = []; argSS = []; argSS2 = []; readSS2 = 0;

if ni>=2,
   if ~isempty(degree),
      if ~isa(degree,'double') | ~isreal(degree) | ...
            ndims(degree)~=2,
         error('Invalid 2nd argument; must be degree.');
      end;
      degree = degree(:); ldeg = length(degree);
      if ldeg>2 | ~all(floor(degree)==degree),
         error('Invalid 2nd argument; must be degree.');
      end;
      if ldeg==1,
         if degree<0 & isfinite(degree),
            error('Invalid 2nd argument; must be degree.');
         end;
         Ad = degree; Atd = -degree;
      else  % ldeg==2
         Atd = degree(1); Ad = degree(2);
      end;
   end;
end;

if ni>=3,
   if ~isempty(T),
      if isa(T,'double') & length(T)==1 & isreal(T) & T>=0,
         Ah = T;
      else
         error('Invalid 3rd argument; must be sampling period.');
      end;
   end;
end;

if ni>=4,
   if ~isempty(tau),
      if isa(tau,'double') & ndims(tau)==2 & any(size(tau)==1) & ...
            isreal(tau) & all(tau>=0),
         Atau = sort(tau);
      else
         error('Invalid 4th argument; must be sampling phase.');
      end;
   end;
end;

lv = length(varargin);
for i = 1:lv,
   arg = varargin{i};
   if isa(arg,'double'),
      if ~isempty(arg),
         if ndims(arg)~=2,
            error(['Invalid ',nth(i+4),' argument.']);
         end;
         args = size(arg);
         if ~all(args==As),
            error(['Invalid ',nth(i+4), ...
                  ' argument; matrices not of the same dimensions.']);
         end;
         if isreal(arg) & all(all(arg>0)) & ...
               all(all(floor(arg)==arg)),
            HH = arg;
         else
            error(['Invalid ', nth(i+4), ' argument.']);
         end;
      end;
   elseif isa(arg,'char'),
      if ~readSS2,
         argSS = arg; readSS2 = 1;
      else
         argSS2 = arg;
      end;
   else
      error(['Invalid ', nth(i+4), ' argument.']);
   end;
end;

if isempty(Ah) | ~isfinite(Ah) | Ah==0,
   Ah = 1; SS = 'ks'; SS2 = 'kd';
else
   SS = 'bs'; SS2 = 'bd';
end;

if Atd==-inf, Atd = Atdeg;
end;
if Ad==+inf, Ad = Adeg;
end;

if Ad<Atd,
   Ac = zeros(As1,As2,1); Ad = 0; Atd = 0; As3 = 1;
   SS = 'ks'; SS2 = 'kd';
elseif ~isfinite(Adeg),
   Ac = zeros(As1,As2,Ad-Atd+1);
   As3 = Ad-Atd+1;
else
   if Atd>Adeg | Ad<Atdeg,
      Ac = zeros(As1,As2,Ad-Atd+1);
      As3 = Ad-Atd+1;
   else
      Ac = A.p.c; As3 = size(Ac,3);
      Ac = Ac(:,:,As3:-1:1);
      if Atd<Atdeg,
         Ac = cat(3,zeros(As1,As2,Atdeg-Atd),Ac);
         As3 = As3+Atdeg-Atd;
      elseif Atd>Atdeg,
         Ac = Ac(:,:,Atd-Atdeg+1:As3);
         As3 = As3-Atd+Atdeg;
      end;
      if Ad>Adeg,
         Ac = cat(3,Ac,zeros(As1,As2,Ad-Adeg));
         As3 = As3+Ad-Adeg;
      elseif Ad<Adeg,
         Ac = Ac(:,:,1:As3-Adeg+Ad);
         As3 = As3-Adeg+Ad;
      end;
   end;
end;

if isempty(Ah) | ~isfinite(Ah) | Ah==0,
   Ah = 1; SS = 'ks'; SS2 = 'kd';
end;

if ~isempty(argSS),
   SS = argSS;
end;
if ~isempty(argSS2),
   SS2 = argSS2;
end;

ltau = length(Atau); onesltau = ones(1,ltau);
D = kron(Atd*Ah:Ah:Ad*Ah,onesltau) + repmat(Atau,1,As3);
if isreal(A),
   for i = 1:As1,
      for j = 1:As2,
         if isempty(HH), G = figure;
         else G = figure(HH(i,j));
         end;
         if no, H(i,j) = G;
         end;
         E = kron(reshape(Ac(i,j,1:As3),1,As3),onesltau);
         eval('plot(D,E,SS);', 'error(lasterr);');
      end;
   end;
else
   Acre = real(Ac); Acim = imag(Ac);
   for i = 1:As1,
      for j = 1:As2,
         if isempty(HH), G = figure;
         else G = figure(HH(i,j));
         end;
         if no, H(i,j) = G;
         end;
         Ere = kron(reshape(Acre(i,j,1:As3),1,As3),onesltau);
         Eim = kron(reshape(Acim(i,j,1:As3),1,As3),onesltau);
         eval('plot(D,Ere,SS,D,Eim,SS2);', 'error(lasterr);');
      end;
   end;
end;

%end .. @tsp/picture
         