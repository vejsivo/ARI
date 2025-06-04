function H = picture(A,degree,T,tau,varargin)
%PICTURE    Picture fraction
%
% For scalar fraction in variable 'z','z^-1','d','q',
% the command  PICTURE(A)  pictures A, i.e. plots the
% coefficients of its Laurent series A{K} versus K,
% in powers of z^-1 (to represent time signals in
% two-sided Z-transform).
% The command  PICTURE(A,DEGREE,T,TAU) plots the
% coefficients A{K} up to  K = DEGREE  versus time
%  K*T + TAU . If DEGREE is vector  [DEG1,DEG2] then
% the range of K is DEG1:DEG2. Arguments DEGREE,T,TAU
% may be omitted or given by [] . Default for DEGREE is
%  [ TDEG(Q), MAX(DEG(A.N),DEG(A.D)) ]  where Q is
% two-sided polynomial of Laurent coefficients of A,
% considered in powers of z^-1. Default for T is
%  A.h or 1, for TAU is 0.
%
% For strictly proper scalar fraction in variable 's',
% 'p', the plotted coefficients A{K} correspond to the
% inverse Laplace transform of A, sampled with period T
% and phase TAU, i.e. in sampling points  K*T + TAU,
%  K = 0,1,2...DEGREE. In this case, argument T must
% be given, argument TAU may be [] or scalar or vector.
% If A is not strictly proper, its polynomial part is
% also plotted, in the same figure.
%
% For matrix fraction A, the command pictures
% all entries  A(I,J)  as said above, each entry
% in a separate figure.
%
% Standardly, new figures are created and their
% handles H(I,J), natural numbers, may be returned
% in the output argument  H = PICTURE(A) . 
% However, when an optional input argument H,
% matrix of natural numbers, is used in
%  PICTURE(A,T,TAU,DEGREE,H) , the existing
% figures may be used.
%
% By optional input argument S, character string,
% various line types, plot symbols and colors may
% be specified, for all entries of A uniformly,
% as in PLOT. For discrete-time fractions with
% sampling period, the default is 'bs' (blue square),
% without sampling period is 'ks' (black square).
% For continuous-time fractions, the default is '-'
% (solid line).
%
% If A is complex then the real and the imaginary 
% parts are pictured in the same figure. The default
% for discrete fractions: the real part by squares,
% the imaginary part by diamonds. The default for 
% continuous time fractions: the real part by solid
% line, the imaginary part by dotted one. Instead of
% optional input argument S, two arguments S1,S2 may
% be used.
%
% See also POL/PICTURE, TSP/PICTURE, PLOT, FIGURE.

%         Author:  J.Jezek, 06-Jan-2001
%         Copyright(c) 2001 by Polyx, Ltd.
%         $ Revision $  $ Date 25-Jul-2002 $
%                       $ Date 14-Oct-2002 $

ni = nargin; no = nargout;
if ni<1,
   error('Not enough input arguments.');
end;
Acl = class(A);
if strcmp(Acl,'frac'),
   error('Invalid 1st argument.');
end;
if ~isa(A,'frac'),
   error('Invalid 1st argument.');
end;

if isempty(A),
   if no, H = zeros(size(A));
   end;
   return;
end;

As = size(A); As1 = As(1); As2 = As(2);
Ah = A.h; Atau = 0;
Adeg = deg(A.num);
if isfinite(Adeg),
   Adeg = max(Adeg, deg(A.den));
end;
Ad = Adeg; Atd = []; dirac = 0;
HH = []; argSS = []; argSS2 =[]; readSS2 = 0;

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
         Ad = degree;
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
ltau = length(Atau); stau = size(Atau);

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

Av = A.v;
if isempty(Av) | strcmp(Av,'z^-1') | strcmp(Av,'d') | ...
      strcmp(Av,'z') | strcmp(Av,'q'),  % discrete time
   
   if isempty(Ah) | ~isfinite(Ah) | Ah==0,
      Ah = 1; SS = 'ks'; SS2 = 'kd';
   else
      SS = 'bs'; SS2 = 'bd';
   end;

   if Ad==+inf, Ad = Adeg;
   end;
   
   if Ad~=-inf,
      if Ad>=0,
         Q = laurent(A,Ad);
         if strcmp(Q.v,'d'), Q.v = 'zi';
         end;
         Q = tsp(Q);
      else
         Q = pol(laurent(A,0));
         Q = shift(Q,Ad);
         Q = tsp(shift(Q,-Ad));
      end;
      Qd = -tdeg(Q);
      Ac = Q.p.c; As3 = size(Ac,3);
      Ac = Ac(:,:,As3:-1:1);
      if isfinite(Qd) & Qd<Ad,
         Ac = cat(3,Ac,zeros(As1,As2,Ad-Qd));
         As3 = As3+Ad-Qd;
      end;
      Atdeg = -deg(Q);
   else
      Atdeg = +inf;
   end;
     
   if isempty(Atd) | Atd==-inf, Atd = Atdeg;
   end;
   if Atd>Ad | Atd==+inf | Ad==-inf,
      Ac = zeros(As1,As2,1); Ad = 0; Atd = 0; As3 = 1;
      SS = 'ks'; SS2 = 'kd'; ltau = 1;
   elseif Atdeg==+inf | Adeg==-inf,
      Ac = zeros(As1,As2,Ad-Atd+1);
      As3 = Ad-Atd+1;
   else
      As3 = Ad-Atdeg+1;
      if Atd<Atdeg,
         Ac = cat(3,zeros(As1,As2,Atdeg-Atd),Ac);
         As3 = As3+Atdeg-Atd;
      elseif Atd>Atdeg,
         Ac = Ac(:,:,Atd-Atdeg+1:As3);
         As3 = As3-Atd+Atdeg;
      end;
   end;
   
   if ltau~=1,
      Acc = Ac;
      Ac = cell(stau);
      for k = 1:ltau,
         Ac{k} = Acc;
      end;
   end;
   
else    % continuous time
   
   if Ah==0, Ah = 1;
   end;
   SS = 'b-'; SS2 = 'b:'; SS3 = 'ks';  SS4 = 'kd';
   if Ad==+inf, Ad = Adeg;
   elseif Ad<0, Ad = 0;
   end;
   if isempty(Atd) | Atd<0, Atd = 0;
   end;
   
   if Atd>Ad,
      Ac = zeros(As1,As2,1); Ad = 0; Atd = 0; As3 = 1;
      SS = 'ks'; SS2 = 'kd'; ltau = 1;
   else
      [R,Q] = laplace(A,Ad,Ah,Atau);
      
      if ltau==1,
         Rd = deg(R);
         if ~isfinite(Rd),
            Ac = zeros(As1,As2,Ad-Atd+1);
         else
            Ac = R.c;
            if Atd>0,
               Ac = Ac(:,:,Atd+1:Ad+1);
            end;
         end;
         if Atd==0 & Q~=0,
            QQ = pol(Q); Qc = QQ.c; dirac = 1;
            Qs3 = size(Qc,3);
         end;
         As3 = size(Ac,3);
         
      else
         Ac = cell(stau); Qc = Ac;
         for k = 1:ltau,
            Rd = deg(R{k});
            if ~isfinite(Rd),
               Ac{k} = zeros(As1,As2,Ad-Atd+1);
            else
               Ac{k} = R{k}.c;
               if Atd>0,
                  Ac{k} = Ac{k}(:,:,Atd+1,Ad+1);
               end;
            end;
            if Atd==0 & Q{k}~=0,
               QQ = pol(Q{k}); Qc{k} = QQ.c; dirac = 1;
               Qs3 = size(Qc{1},3);
            end;
         end;
         As3 = size(Ac{1},3);
      end;  
   end;
  
end;   

if ~isempty(argSS),
   SS = argSS;
end;
if ~isempty(argSS2),
   SS2 = argSS2;
end;

onesltau = ones(1,ltau);
D = kron(Atd*Ah:Ah:Ad*Ah,onesltau) + repmat(Atau,1,As3);
if ~dirac,
   if isreal(A),
      E = zeros(ltau,As3);
      for i = 1:As1,
         for j = 1:As2,
            if isempty(HH), G = figure;
            else G = figure(HH(i,j));
            end;
            if no, H(i,j) = G;
            end;
            if ltau==1,
               E = reshape(Ac(i,j,1:As3),1,As3);
               eval('plot(D,E,SS);', 'error(lasterr);');
            else
               for k = 1:ltau,
                  E(k,:) = reshape(Ac{k}(i,j,1:As3),1,As3);
               end;
               eval('plot(D,E(:),SS);', 'error(lasterr);');
            end;
         end;
      end;
   else    % ~isreal(A)
      
      Ere = zeros(ltau,As3); Eim = Ere;
      for i = 1:As1,
         for j = 1:As2,
            if isempty(HH), G = figure;
            else G = figure(HH(i,j));
            end;
            if no, H(i,j) = G;
            end;
            if ltau==1,
               Ere = reshape(real(Ac(i,j,1:As3)),1,As3);
               Eim = reshape(imag(Ac(i,j,1:As3)),1,As3);
               eval('plot(D,Ere,SS,D,Eim,SS2);', 'error(lasterr);');
            else
               for k = 1:ltau,
                  Ere(k,:) = reshape(real(Ac{k}(i,j,1:As3)),1,As3);
                  Eim(k,:) = reshape(imag(Ac{k}(i,j,1:As3)),1,As3);
               end;
            end;
            eval('plot(D,Ere(:),SS,D,Eim(:),SS2);', 'error(lasterr);');
         end;   
      end;
   end;
   
else    % dirac
   N = 0 : Ah : Ah*deg(QQ);
   if isreal(A),
      E = zeros(ltau,As3); M = E;
      for i = 1:As1,
         for j = 1:As2,
            if isempty(HH), G = figure;
            else G = figure(HH(i,j));
            end;
            if no, H(i,j) = G;
            end;
            if ltau==1,
               E = reshape(Ac(i,j,1:As3),1,As3);
               M = reshape(Qc(i,j,1:Qs3),1,Qs3);
               eval('plot(D,E,SS,N,M,SS3);', 'error(lasterr);');
            else
               for k = 1:ltau,
                  E(k,:) = reshape(Ac{k}(i,j,1:As3),1,As3);
                  M(k,:) = reshape(Qc{k}(i,j,1:Qs3),1,Qs3);
               end;
               eval('plot(D,E(:),SS,N,M(:),SS3);', 'error(lasterr);');
            end;
         end;
      end;
      
   else   % dirac, ~isreal(A)
      Ere = zeros(ltau,As3); Eim = Ere;
      Mre = Ere; Mim = Ere;
      for i = 1:As1,
         for j = 1:As2,
            if isempty(HH), G = figure;
            else G = figure(HH(i,j));
            end;
            if no, H(i,j) = G;
            end;
            if ltau==1,
               Ere = reshape(real(Ac(i,j,1:As3)),1,As3);
               Eim = reshape(imag(Ac(i,j,1:As3)),1,As3);
               Mre = reshape(real(Qc(i,j,1:Qs3)),1,Qs3);
               Mim = reshape(imag(Qc(i,j,1:Qs3)),1,Qs3);
               eval('plot(D,Ere,SS,D,Eim,SS2,N,Mre,SS3,N,Mim,SS4);',...
                 'error(lasterr);');
            else
               for k = 1:ltau,
                  Ere(k,:) = reshape(real(Ac{k}(i.j,1:As3)),1,As3);
                  Eim(k,:) = reshape(imag(Ac{k}(i,j,1:As3)),1,As3);
                  Mre(k,:) = reshape(real(Qc{k}(i,j,1:As3)),1,As3);
                  Mim(k,:) = reshape(imag(Qc{k}(i,j,1:As3)),1,As3);
               end;
               eval('plot(D,Ere(:),SS,D,Eim(:),SS2,N,Mre(:),SS3,N,Mim(:),SS$);',...
                  'error(lasterr);');
            end;   
         end;
      end;
   end;
      
end;

%end .. @frac/picture
