function P = subsref(R,S)
%SUBSREF    Subscripted reference for scalar-den fraction
%
% For scalar-denominator fraction R,
%   R.num  returns the numerator 
%   R.den          the denominator
%   R.s            the size
%   R.d            the degree, same as  DEG(R)
%   R.v            the variable symbol
%   R.h            the sampling period
%   R.u            the user data
%   R.c            the coprime flag
%   R.r            the reduced flag
%   R.p            the proper  flag
%   R.tc           the coprime tolerance
%   R.tp           the proper tolerance
%   R.version      the version, always 3.0
%   R(i,j)         the scalar fraction
%                    equal to (i,j)-th entry of R
%   R(i1:i2,:)     similarly
%   R(X) where X(var) is scalar polynomial, two-sided pol
%        or fraction, returns the composition  R(X(var)).
%   R{k} returns   the k-th coefficient in Laurent series of R
%                  (at k-th power of F.v)
%
% For the numerator and denominator, the subscripting in reference
% may be further elaborated, as for polynomials, e.g.  R.num(i,j)
% or  R.den{k}  or  R{k}(i,j) 

%     Author:   J.Jezek  24-Jan-2000
%     Copyright(c) 2000 by Polyx, Ltd.
%     $ Revision $  $ Date 09-May-2000 $
%                   $ Date 22-May-2000 $
%                   $ Date 20-Oct-2000 $
%                   $ Date 06-Feb-2002 $
%                   $ Date 17-Jul-2002 $
%                   $ Date 06-Oct-2002 $
%                   $ Date 14-Oct-2002 $

global PGLOBAL;

done = 1;
St1 = S(1).type; Ssubs = S(1).subs;
if strcmp(St1,'.') & strcmp(Ssubs,'frac'),
   S = S(2:end);
   St1 = S(1).type; Ssubs = S(1).subs;
end;
if strcmp(St1,'.'),
   if strcmp(Ssubs,'version'), P = R.frac.version;
   elseif strcmp(Ssubs,'den'), P = R.frac.den;            
   else   
      switch lower(Ssubs(1)),
      case 'n', P = R.frac.n;
      case 's', P = R.frac.s;
      case 'd', P = deg(R);
      case 'v', P = R.frac.v;
      case 'h', P = R.frac.h;
      case 'u', P = R.frac.u;
      case 'c', P = R.frac.c;
      case 'r', P = R.frac.r;
      case 'p', P = R.frac.p;
      case 't',
         if strcmp(Ssubs,'tc'), P = R.frac.tc;
         elseif strcmp(Ssubs,'tp'), P = R.frac.tp;
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
      if isa(X,'pol') | isa(X,'tsp') | isa(X,'frac'),
         if any(size(X)~=1),
            error('Invalid argument; must be scalar polynomial or fraction.');
         end;
         n = max(deg(R.frac.num), deg(R.frac.den));
         if ~isempty(n) & isfinite(n),
            if isa(X,'pol') | isa(X,'tsp'),
               NN = zeros(size(R)); DD = 0;
               for k = n:-1:0,
                  NN = NN*X + R.frac.num{k};
                  DD = DD*X + R.frac.den{k};
               end;
            else
               NN = R.frac.num{0}; H = 1;
               DD = R.frac.den{0};
               for k = 1:n,
                  H = H*X.frac.num;
                  NN = NN*X.frac.den + R.frac.num{k}*H;
                  DD = DD*X.frac.den + R.frac.den{k}*H;
               end;
            end;
         else
            NN = R.frac.num; DD = R.frac.den;
         end;
         P = sdf(NN,DD);
      else
         i = Ssubs{1};
         NN = 0; DD = 0;
         eval('NN = R.frac.num(i); DD = R.frac.den;','error(peel(lasterr));');
         P = sdf(NN,DD);
         if strcmp(R.frac.p,'prop'), props(P,'prop',R.frac.tp);
         end;
         if strcmp(R.frac.r,'red'), props(P,'red');
         end;
      end;
   elseif lS==2,
      i = Ssubs{1}; j = Ssubs{2};
      NN = 0; DD = 0;
      eval('NN = R.frac.num(i,j); DD = R.frac.den;','error(peel(lasterr));');
      P = sdf(NN,DD);
      if strcmp(R.frac.p,'prop'), props(P,'prop',R.frac.tp);
      end;
      if strcmp(R.frac.r,'red'), props(P,'red');
      end;
   else
      error('More than two subscripts.');
   end;
   if strcmp(PGLOBAL.COPRIME,'cop'), P = coprime(P);
   end;
   if strcmp(PGLOBAL.REDUCE,'red'), P = reduce(P);
   else P = smreduce(P);
   end;
   if strcmp(PGLOBAL.DEFRACT,'defr'), P = defract(P);
   end;
elseif strcmp(St1,'{}'),
   lS = length(Ssubs);
   if lS==1,
      k = Ssubs{1};
      if ~isa(k,'double') | ~isreal(k),
         error('Invalid subscript in {} ; must be integer.');
      end;
      if strcmp(R.frac.v,'s') | strcmp(R.frac.v,'p') | ...
            strcmp(R.frac.v,'q') | strcmp(R.frac.v,'z'),
         R.frac.v = 'z';
         k = -k; 
      elseif strcmp(R.frac.v,'d'),
         R.frac.v = 'z^-1';
      end;      
      maxk = max(max(k));
      if ~isfinite(maxk),
         error('Subscript is not finite.');
      end;      
      maxk = max(maxk,0);
      H = tsp(laurent(R,maxk));
      P = H{-k};
   else
      error('More than one subscript in {}.');
   end; 
   if length(S)>=2,
      St2 = S(2).type;
      if strcmp(St2,'()'),
         Ssubs = S(2).subs; lS = length(Ssubs);
         Rs1 = R.frac.s(1); Rs2 = R.frac.s(2);Rs12 = Rs1*Rs2;
         [Ps1,Ps2] = size(P); Ps12 = Ps1*Ps2;
         if lS==1,
            i = Ssubs{1};
            P = reshape(P,Rs12,Ps12/Rs12);
            eval('P = P(i,:);','error(lasterr);');
            [Ps1,Ps2] = size(P);
            P = reshape(P,1,Ps1*Ps2);
         elseif lS==2,
            i = Ssubs{1}; j = Ssubs{2};
            P = reshape(P,Rs1,Rs2,Ps2/Rs2);
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
      
%end .. @sdf/subsref
