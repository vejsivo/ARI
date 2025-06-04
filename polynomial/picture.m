function H = picture(A,varargin)
%PICTURE    Picture constant or cell vector
%
% For constant A (i.e.a standard Matlab matrix),
% the command  PICTURE(A)  pictures A as POL(A).
% See POL/PICTURE, also for possible further
% arguments.
%
% Argument A may also be a cell vector whose cells
% contain discrete-time polynomials, tsp or fractions.
% All cells correspond to the same continuous-time
% object and to the same sampling period T, each cell
% for different sampling phase TAU. The command
%  PICTURE(A,DEGREE,T,TAU) where TAU is a vector of
% sampling phases, pictures A, all TAUs in the same
% figure.
%
% For more details and possible further arguments,
% see also POL/PICTURE, TSP/PICTURE, FRAC/PICTURE.

%        Author :  J.Jezek, 16-Jan-2001
%        Copyright(c) 2001 by Polyx, Ltd.
%        $ Revision $  $ Date 20-Jan-2001 $

ni = nargin; no = nargout;
if ni<1,
   error('Not enough input arguments.');
end;
lv = length(varargin);

if isa(A,'double'),
   if ndims(A)~=2,
      error('Invalid 1st argument.');
   end;
   eval('picture(pol(A),varargin{1:lv});', ...
      'error(peel(lasterr));');
elseif isa(A,'cell'),
   if ndims(A)~=2,
      error('Invalid 1st argument.');
   end;
   cs = size(A); Al = length(A);
   if ~any(cs==1) | Al<1,
      error('Invalid 1st argument.');
   end;
   
   A1 = A{1}; 
   if isa(A1,'pol'),
      Adeg = deg(A1); Atdeg = 0;
   elseif isa(A1,'tsp'),
      Adeg = -tdeg(A1); Atdeg = -deg(A1);
   elseif isa(A1,'frac'),
      Adeg = deg(A1.n);
      if isfinite(Adeg),
         Adeg = max(Adeg,deg(A1.d));
      end;
      Atdeg = [];
   else
      error('Invalid 1st argument.');
   end;
      
   Av = A1.v;
   if strcmp(Av,'s') | strcmp(Av,'p'),
      error('Invalid 1st argument; must be discrete time.');
   end;
   Ah = A1.h; Acl = class(A1);
   As = size(A1); As1 = As(1); As2 = As(2);
   if As1==0 | As2==0,
      if no, H = zeros(As);
      end;
      return;
   end;
   
   for k = 2:Al,
      Ak = A{k};
      if ~strcmp(class(Ak),Acl),
         error('Invalid 1st argument; object classes not the same.');
      end;
      if ~all(size(Ak)==As),
         error('Invalid 1st argument; matrix dimensions not the same.');
      end;
      
      Akv = Ak.v;
      if ~isempty(Akv),
         if strcmp(Akv,'s') | strcmp(Akv,'p'),
            error('Invalid first argument; must be dicrete time.');
         end;
         if ~isempty(Av),
            if ~strcmp(Akv,Av),
               error('Invalid 1st argument; variable symbols not the same.');
            end;
         else
            Av = Akv;
         end;
      end;
      
      Akh = Ak.h;
      if ~isempty(Akh) & isfinite(Akh),
         if ~isempty(Ah) & isfinite(Ah),
            if Akh~=Ah,
               error('Invalid 1st argument; sampling periods not the same.');
            end;
         else
            Ah = Akh;
         end;
      end;
      
      if strcmp(Acl,'pol'),
         Akdeg = deg(Ak); Aktdeg = 0;
      elseif strcmp(Acl,'tsp'),
         Akdeg = -tdeg(Ak); Aktdeg = -deg(Ak);
      else   % 'frac'
         Akdeg = deg(Ak.n); Aktdeg = [];
         if isfinite(Akdeg),
            Akdeg = max(Akdeg,deg(Ak.d));
         end;
      end;
      Adeg = max(Adeg,Akdeg);
      Atdeg = min(Atdeg,Aktdeg);
   end;
   
   Ad = Adeg; Atd = Atdeg; Atau = 0;
   HH = []; argSS = []; argSS2 = []; readSS2 = 0;
   
   if lv>=1,
      degree = varargin{1};
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

   if lv>=2,
      T = varargin{2};
      if ~isempty(T),
         if isa(T,'double') & length(T)==1 & isreal(T) & T>=0,
            Ah = T;
         else
            error('Invalid 3rd argument; must be sampling period.');
         end;
      end;
   end;
   
   if lv>=3,
      tau = varargin{3};
      if ~isempty(tau),
         if isa(tau,'double') & ndims(tau)==2 & any(size(tau)==1) & ...
               isreal(tau) & all(tau>=0),
            Atau = sort(tau);
         else
            error('Invalid 4th argument; must be sampling phase.');
         end;
      end;
   end;
   if length(Atau)~=Al,
      error('Invalid length of sampling phase vector.');
   end;

   for i = 4:lv,
      arg = varargin{i};
      if isa(arg,'double'),
         if ~isempty(arg),
            if ndims(arg)~=2,
               error(['Invalid ',nth(i+1),' argument.']);
            end;
            args = size(arg);
            if ~all(args==As),
               error(['Invalid ',nth(i+1), ...
                     ' argument; matrices not of the same dimensions.']);
            end;
            if isreal(arg) & all(all(arg>0)) & ...
                  all(all(floor(arg)==arg)),
               HH = arg;
            else
               error(['Invalid ', nth(i+1), ' argument.']);
            end;
         end;
      elseif isa(arg,'char'),
         if ~readSS2,
            argSS = arg; readSS2 = 1;
         else
            argSS2 = arg;
         end;
      else
         error(['Invalid ', nth(i+1), ' argument.']);
      end;
   end;
   
   if isempty(Ah) | ~isfinite(Ah) | Ah==0,
      Ah = 1; SS = 'ks'; SS2 = 'kd';
   else
      SS = 'bs'; SS2 = 'bd';
   end;
   
   Ac = cell(cs);
   if Ad==+inf, Ad = Adeg;
   end;
   
   if strcmp(Acl,'pol'),
      if Atd<0, Atd = 0;
      end;
      if Atd>Ad,
         Al = 1; Atau = 0; Ac{1} = zeros(As1,As2,1);
         Ad = 0; Atd = 0; As3 = 1;
         SS = 'ks'; SS2 = 'kd';
      else
         for k = 1:Al,
            Ak = A{k};
            Akdeg = deg(Ak);
            Akc = Ak.c;
            if Ad<=Akdeg,
               Ac{k} = Akc(:,:,Atd+1:Ad+1);
            elseif Atd<=Akdeg,
               Ac{k} = cat(3,Akc(:,:,Atd+1:Akdeg+1), ...
                  zeros(As1,As2,Ad-Akdeg));
            else
               Ac{k} = zeros(As1,As2,Ad-Atd+1);
            end;
         end;
         As3 = Ad-Atd+1;
      end;
      
   elseif strcmp(Acl,'tsp'),
      if Atd==-inf, Atd = Atdeg;
      end;
      if Atd>Ad,
         Al = 1; Atau = 0; Ac{1} = zeros(As1,As2,1);
         Ad = 0; Atd = 0; As3 = 1;
         SS = 'ks'; SS2 = 'kd';
      else
         for k = 1:Al,
            Ak = A{k};
            Akdeg = -tdeg(Ak); Aktdeg = -deg(Ak);
            if Atd>Akdeg | Ad<Aktdeg,
               Ac{k} = zeros(As1,As2,Ad-Atd+1);
            else
               Akc = Ak.p.c; Aks3 = size(Akc,3);
               Ac{k} = Akc(:,:,Aks3:-1:1);
               if Atd<Aktdeg,
                  Ac{k}=cat(3,zeros(As1,As2,Aktdeg-Atd),Ac{k});
                  Aks3 = Aks3+Aktdeg-Atd;
               elseif Atd>Aktdeg,
                  Ac{k} = Ac{k}(:,:,Atd-Aktdeg+1:Aks3);
                  Aks3 = Aks3-Atd+Aktdeg;
               end;
               if Ad>Akdeg,
                  Ac{k} = cat(3,Ac{k},zeros(As1,As2,Ad-Akdeg));
                  Aks3 = Aks3+Ad-Akdeg;
               elseif Ad<Akdeg,
                  Ac{k} = Ac{k}(:,:,1:Aks3-Akdeg+Ad);
                  Aks3 = Aks3-Akdeg+Ad;
               end;
            end;
         end;
         As3 = Ad-Atd+1;
      end;
      
   else    %  'frac'
      if Ad==-inf, Atd = +inf;
      else
         Qd = zeros(cs); Qtd = Qd;
         for k = 1:Al,
            Ak = A{k};
            if Ad>=0,
               Qk = laurent(Ak,Ad);
               if strcmp(Qk.v,'d'), Qk.v = 'zi';
               end;               
               Qk = tsp(Qk);
            else
               Qk = pol(laurent(Ak,0));
               Qk = shift(Qk,Ad); 
               Qk = tsp(shift(Qk,-Ad));
            end;
            Qd(k) = -tdeg(Qk); Qtd(k) = -deg(Qk);
            Akc = Qk.p.c; Aks3 = size(Akc,3);
            Ac{k} = Akc(:,:,Aks3:-1:1);
         end;
         Atdeg = min(Qtd);
         if isempty(Atd) | Atd==-inf, Atd = Atdeg;
         end;
      end;
      
      if Atd>Ad,
         Al = 1; Atau = 0; Ac{1} = zeros(As1,As2,1);
         Ad = 0; Atd = 0; As3 = 1;
         SS = 'ks'; SS2 = 'kd';
      else
         for k = 1:Al,
            Qdk = Qd(k); Qtdk = Qtd(k);
            if Qtdk>Qdk,
               Ac{k} = zeros(As1,As2,Ad-Atd+1); 
            else
               As3k = Qdk-Qtdk+1;
               if Qtdk>Atd,
                  Ac{k} = cat(3,zeros(As1,As2,Qtdk-Atd),Ac{k});
                  As3k = As3k+Qtdk-Atd;
               elseif Qtdk<Atd,
                  Ac{k} = Ac{k}(:,:,Atd-Qtdk+1:As3k);
                  As3k = As3k-Atd+Qtdk;
               end;
               if Qdk<Ad,
                  Ac{k} = cat(3,Ac{k},zeros(As1,As2,Ad-Qdk));
               elseif Qdk>Ad,
                  Ac{k} = Ac{k}(:,:,1:As3k-Qdk+Ad);
               end;
            end;
         end;
         As3 = Ad-Atd+1;
      end;
      
   end;    %   pol, tsp, frac
   
   if ~isempty(argSS),
      SS = argSS;
   end;
   if ~isempty(argSS2),
      SS2 = argSS2;
   end;
   
   onesAl = ones(1,Al);
   D = kron(Atd*Ah:Ah:Ad*Ah,onesAl) + repmat(Atau,1,As3);
   if isreal(A1),
      E = zeros(Al,As3);
      for i = 1:As1,
         for j = 1:As2,
            if isempty(HH), G = figure;
            else G = figure(HH(i,j));
            end;
            if no, H(i,j) = G;
            end;
            for k = 1:Al,
               E(k,:) = reshape(Ac{k}(i,j,1:As3),1,As3);
            end;
            eval('plot(D,E(:),SS);', 'error(lasterr);');
         end;
      end;
   else   % ~isreal(A)
      Ere = zeros(Al,As3); Eim = Ere;
      for i = 1:As1,
         for j = 1:As2,
            if isempty(HH), G = figure;
            else G = figure(HH(i,j));
            end;
            if no, H(i,j) = G;
            end;
            for k = 1:Al,
               Ere(k,:) = reshape(real(Ac{k}(i,j,1:As3)),1,As3);
               Eim(k,:) = reshape(imag(Ac{k}(i,j,1:As3)),1,As3);
            end;
            eval('plot(D,Ere(:),SS,D,Eim(:),SS2);', 'error(lasterr);');
         end;
      end;
      
   end;
   
else    %  ~(isa(A,'double') | isa(A,'cell'))
   error('Invalid 1st argument.');
end;

%end .. picture

