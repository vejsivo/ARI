function P = subsref(F,S)
%SUBSREF    Subscripted reference for left-den fraction
%
% For left-denominator fraction F,
%   F.num returns the numerator 
%   F.den         the denominator
%   F.s           the size
%   F.d           the degree, same as  DEG(F)
%   F.v           the variable symbol
%   F.h           the sampling period
%   F.u           the user data
%   F.version     the version, always 3.0
%   F.c           the coprime flag
%   F.r           the reduced flag
%   F.p           the proper  flag
%   F.tc          the coprime tolerance
%   F.tp          the proper  tolerance
%   F(i,j)        the right fraction equal to (i,j)-th entry of 
%   F(i,j)        the right fraction equal to (i,j)-th entry of F
%   F(i1:i2,:)    similarly
%   F(X) where X(var) is polynomial, two-sided pol or fraction,
%        returns  the composition  F(X(var)).
%                 If F is scalar, X must be square matrix.
%                 If not, X must be scalar.
%   F{k} returns  the k-th coefficient in Laurent series of F
%                 (at k-th power of F.v)
%
% For the numerator and denominator, the subscripting in reference
% may be further elaborated, as for polynomials, e.g.  F.num(i,j)
% or  F.den{k}  or  F{k}(i,j) 

%     Author:   J.Jezek  02-Dec-1999
%     Copyright(c) 1999 by Polyx, Ltd.
%     $ Revision $  $ Date 09-May-2000 $
%                   $ Date 06-Feb-2002 $
%                   $ Date 17-Jul-2002 $
%                   $ Date 06-Oct-2002 $
%                   $ Date 14-Oct-2002 $
%                   $ Date 05-Nov-2002 $
%                   $ Date 28-Feb-2003 $

global PGLOBAL;

done = 1;
St1 = S(1).type; Ssubs = S(1).subs;
if strcmp(St1,'.') & strcmp(Ssubs,'frac'),
   S = S(2:end);
   St1 = S(1).type; Ssubs = S(1).subs;
end;
if strcmp(St1,'.'),
   if strcmp(Ssubs,'version'), P = F.frac.version;
   elseif strcmp(Ssubs,'den'), P = F.frac.den;
   else   
      switch lower(Ssubs(1)),
      case 'n', P = F.frac.num;
      case 's', P = F.frac.s;
      case 'd', P = deg(F);
      case 'v', P = F.frac.v;
      case 'h', P = F.frac.h;
      case 'u', P = F.frac.u;
      case 'c', P = F.frac.c;
      case 'r', P = F.frac.r;
      case 'p', P = F.frac.p;
      case 't',
         if strcmp(Ssubs,'tc'), P = F.frac.tc;
         elseif strcmp(Ssubs,'tp'), P = F.frac.tp;
         else error('Invalid field in subscripted reference.');
         end;         
      otherwise
         error('Invalid field in subscripted reference.');
      end;
   end;
elseif strcmp(St1,'()'),
   lS = length(Ssubs);
   if lS==1,
      X = Ssubs{1};  % pol, tsp or frac in place of subscript
      if isa(X,'tsp'), X = sdf(X);
      end;
      if isa(X,'pol') | isa(X,'frac'),
         if any(size(F)~=1),
            if any(size(X)~=1),
               error('Invalid argument; must be scalar polynomial or fraction.');
            end;
         else
            if size(X,1)~=size(X,2),
               error('Invalid argument; must be square matrix pol or frac.');
            end;
         end;
         n = max(deg(F.frac.num), deg(F.frac.den));
         if ~isempty(n) & isfinite(n),
            if isa(X,'pol') | isa(X,'tsp') | any(size(X)~=1),
               e = eye(size(X));
               NN = zeros(size(F.frac.num));
               DD = zeros(size(F.frac.den));
               for k = n:-1:0,
                  NN = NN*X + F.frac.num{k}*e;
                  DD = DD*X + F.frac.den{k}*e;
               end;
               P = ldf(DD,NN);
            else
               NN = F.frac.num{0}; H = 1;
               DD = F.frac.den{0};
               for k = 1:n,
                  H = H*X.frac.num;
                  NN = NN*X.frac.den + F.frac.num{k}*H;
                  DD = DD*X.frac.den + F.frac.den{k}*H;
               end;
               P = ldf(DD,NN);
            end;
         else
            P = F;
         end;
         if strcmp(PGLOBAL.COPRIME,'cop'), P = coprime(P);
         end;
         if strcmp(PGLOBAL.REDUCE,'red'), P = reduce(P);
         else P = smreduce(P);
         end;
         if strcmp(PGLOBAL.DEFRACT,'defr'), P = defract(P);
         end;
      else
         FF = reshape(F,size(F,1)*size(F,2),1);
         S(1).subs{1} = Ssubs{1};
         S(1).subs{2} = 1;
         eval('P = subsref(FF,S);','error(peel(lasterr));');
         return;
      end;
   elseif lS==2,
      Fp = F.frac.p; j = Ssubs{2};
      eval('F.frac.num = F.frac.num(:,j);','error(peel(lasterr));');
      i = Ssubs{1}; n = size(F.frac.den,1);
      if isa(i,'char'),
         if ~strcmp(i,':'),
            error('Invalid 1st subscript.');
         end;
         e = eye(n);
      else                     %  isa(i,'cell')
         li = length(i);
         e = zeros(li,n);
         for k = 1:li,
            ik=i(k);
            if ~any(ik==(1:n)),
               error('1st subscript invalid or exceeds matrix dimensions.');
            end;
            e(k,ik) = 1;
         end;
      end;
      P = e*ldf(F.frac.den,F.frac.num);
      if strcmp(Fp,'prop'),props(P,'prop',F.frac.tp);
      end;
   else
      error('More than two subscripts.');
   end;
elseif strcmp(St1,'{}'),
   lS = length(Ssubs);
   if lS==1,
      k = Ssubs{1};
      if ~isa(k,'double') | ~isreal(k),
         error('Invalid subscript in {} ; must be integer.');
      end;
      if strcmp(F.frac.v,'s') | strcmp(F.frac.v,'p') | ...
            strcmp(F.frac.v,'q') | strcmp(F.frac.v,'z'),
         F.frac.v = 'z';
         k = -k; 
      elseif strcmp(F.frac.v,'d'),
         F.frac.v = 'z^-1';
      end;
      maxk = max(max(k));
      if ~isfinite(maxk),
         error('Subscript is not finite.');
      end      
      maxk = max(maxk,0);
      H = tsp(laurent(F,maxk));
      P = H{-k};
   else
      error('More than one subscript in {} .');
   end; 
   if length(S)>=2,
      St2 = S(2).type;
      if strcmp(St2,'()'),
         Ssubs = S(2).subs; lS = length(Ssubs);
         Fs1 = F.frac.s(1); Fs2 = F.frac.s(2); Fs12 = Fs1*Fs2;
         [Ps1,Ps2] = size(P); Ps12 = Ps1*Ps2;
         if lS==1,
            i = Ssubs{1};
            P = reshape(P,Fs12,Ps12/Fs12);
            eval('P = P(i,:);','error(lasterr);');
            [Ps1,Ps2] = size(P);
            P = reshape(P,1,Ps1*Ps2);
         elseif lS==2,
            i = Ssubs{1}; j = Ssubs{2};
            P = reshape(P,Fs1,Fs2,Ps2/Fs2);
            eval('P = P(i,j,:);','error(lasterr);');
            [Ps1,Ps2,Ps3] = size(P);
            P = reshape(P,Ps1,Ps2*Ps3);
         else
            error('More than two subscripts.');
         end;
         done = 2;
      end;
   end;
end;
if length(S)>done,
   if (isa(P,'double') | isa(P,'char')) & strcmp(S(done+1).type,'.'),
      error('Invalid combination of subscripted references.');
   end;
   eval('P = subsref(P,S(done+1:end));', ...
      ['if isa(P,''pol'') | isa(P,''tsp'') | isa(P,''frac''),', ...
         'error(peel(lasterr)); else error(peelf(lasterr)); end;'])   
end;
      
%end .. @ldf/subsref
