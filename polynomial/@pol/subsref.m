function P = subsref(P,S)
%SUBSREF   Subscripted reference for polynomial
%
% If P is the polynomial matrix 
%   P0 + P1*s + P2*s^2 + .. + Pd*s^d
% then
%    P{di} 
% returns the block row constant matrix corresponding to 
% the degree indices di, in particular
%    P{0:d} is [P0 P1 P2 .. Pd],
%    P{0} is P0, the constant term of P, and
%    P{k} is Pk, the  k-term of P.
%    P(ri, ci) returns the polynomial submatrix of P of row 
%              indices ri and column indices ci.
%    P(X) where X(var) is polynomial, two-sided pol or fraction,
%         returns the composition  P(X(var)).
%         If P is scalar polynomial, X must be square matrix.
%         If not, X must be scalar.
%
%    P.deg  returns the degree of P.
%    P.size returns the size of P.
%    P.coef returns the three dimensional coefficient matrix of P.
%    P.var  returns the indeterminate string variable of P.
%    P.h    returns the sampling period of P.
%    P.user returns the user data of P.
%
% The above subscripted references may be combined, i.e.
% P(ri,ci).deg returns the degree of a submatrix of P, and
% P{di}(ri,ci) returns the constant coefficient matrix
% corresponding to the degree di of the submatrix P(ri,ci).
%
% See also POL/SUBSASGN.

%       Author(s):  S. Pejchova, D. Henrion, M. Sebek 25-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 05-Nov-1999 09:55:34   $
%       $Revision: 3.0 $  $Date: 04-May-2000  J.Jezek   $ 
%                         $Date: 17-Aug-2000  J.Jezek   $
%                         $Date: 30-Nov-2000  J.Jezek   $
%                         $Date: 24-Jun-2001  J.Jezek   $
%                         $Date: 01-Feb-2002  J.Jezek   $
%                         $Date: 17-Jul-2002  J.Jezek   $
%                         $Date: 03-Nov-2002  J.Jezek   $
%                         $Date: 28-Feb-2003  J.Jezek   $

if nargin<2,
   error('Not enough input arguments.');
end;
if ~isstruct(S),
   error('Invalid 2nd argument.');
end;

P_size = size(P); Pd = deg(P);
S_lg = length(S); S1 = S(1);
St1 = S1.type; Ss1 = S1.subs;
Ss1_lg = length(Ss1);
switch St1,
case '.',
   flstr = lower(Ss1); done = 0;
   if length(flstr)>=3,
      flstr = flstr(1:3);
      if strcmp(flstr,'ver'),
         P = P.version; done = 1;
      elseif strcmp(flstr,'den'),
         P = pol(1); done = 1;
      end;
   end;
   if ~done,
      flstr1 = flstr(1);
      switch flstr1,
      case 'f',
         if S_lg==1,
            error('Invalid field in subscripted reference.');
         end;
         S(1) = []; S_lg = length(S);
         St1 = S(1).type; Ss1 = S(1).subs;
         if ~strcmp(St1,'.'),
            error('Invalid combination in subscripted reference.');
         end;
         flstr = lower(Ss1); flstr1 = flstr(1);
         switch flstr1,
         case 'n',
         case 'd',
            P = pol(1);
         otherwise,
            error('Invalid field in subscripted reference.');
         end;   
      case 'n',
      case 'd', P = P.d;
      case 'v', P = P.v;
      case 's', P = P.s;
      case 'c', P = P.c;
      case 'u', P = P.u;
      case 'h', P = P.h;
      otherwise,
         error('Invalid field in subscripted reference.');
      end;
   end;
   
case '()',
   if Ss1_lg >2,
      error('More than two subscripts.');
   end;
   Sp = Ss1{1};
   if Ss1_lg==1 & (isa(Sp,'pol') | isa(Sp,'tsp') | isa(Sp,'frac')),
      Sp_size = size(Sp);
      if any(P_size~=1),
         if any(Sp_size~=1),
            error('Invalid argument; must be scalar polynomial or fraction.');
         end;
      else
         if Sp_size(1)~=Sp_size(2),
            error('Invalid argument; must be square matrix polynomial or fraction.');
         end;
      end;
      if ~isempty(Pd) & isfinite(Pd),
         if isa(Sp,'pol') | isa(Sp,'tsp') | any(Sp_size~=1),
            e = eye(Sp_size);
            G = zeros(Sp_size);
            for k = Pd:-1:0,
               G = G*Sp + P.c(:,:,k+1)*e;
            end;
            P = G;
         else
            G = P.c(:,:,1); H = 1;
            for k = 1:Pd,
               H = H*Sp.num;
               G = G*Sp.den + P.c(:,:,k+1)*H;
            end;
            P = repmat(Sp,P_size);
            P.num = G; P.den = Sp.den^Pd;
         end;
      end;
         
   else
      if isempty(Pd) | isinf(Pd),
         Pd = 0; P.c = zeros(P_size);
      end;
      if Ss1_lg==2,
         Ppc = P.c; S1.subs{3} = ':';
      elseif isempty(Sp),
         Ppc = zeros(P_size);
      else
         Ppc = reshape(P.c,[P_size(1)*P_size(2),1,Pd+1]);
         S1.subs{2} = [1]; S1.subs{3} = ':';
      end;
      if Ss1_lg>=1,
         i1 = Ss1{1};
         if isnumeric(i1) & ~isempty(i1) & i1<=0,
            error('Index into matrix is negative or zero.');
         end;
      end;
      if Ss1_lg>=2,
         i2 = Ss1{2};
         if isnumeric(i2) & ~isempty(i2) & 12<=0,
            error('Index into matrix is negative or zero.');
         end;
      end;
      Px = 0;
      eval('Px = subsref(Ppc,S1);','error(peelf(lasterr));');
      [ps1,ps2,ps3] = size(Px);
      P.d = ps3-1;
      if ~ps1 | ~ps2,
         P.d = []; Px = zeros(ps1,ps2,0); P.v = '';
      end;
      P.s = [ps1,ps2];
      P.c = Px;
      P.u = [];
      P = pclear(P);
   end;
   
case '{}',
   if Ss1_lg>1,
      error('More than one subscript in {} .');
   end;
   Sp = Ss1{1};
   S1.type = '()';
   S1.subs{1} = ':'; S1.subs{2} = ':'; S1.subs{3} = Sp;
   if isempty(Pd) | isinf(Pd),
      Pd = 0; Ppc = zeros(P_size);
      Sp = zeros(size(Sp));
   else
      Ppc = P.c;
   end;
   if isa(Sp,'double'),
      if any(any(~isfinite(Sp))),
         error('Subscript is not finite.');
      end;
      if any(any(Sp<0)),
%         Pv = P.v;
%         if isempty(Pv) | strcmp(Pv,'z') | strcmp(Pv,'z^-1'),
            Spmin = min(min(Sp));
            Ppc = cat(3,zeros(P_size(1),P_size(2),-Spmin),Ppc);
            Sp = Sp-Spmin;
%         else
%            error('Subscript in {} is negative.');
%         end;
      end;
      
      Spmax = max(max(Sp));
      if Spmax>Pd,
         Ppc = cat(3,Ppc,zeros(P_size(1),P_size(2),Spmax-Pd));
      end;
      S1.subs{3} = Sp+1;
   end;
   eval('Ppc = subsref(Ppc,S1);','error(peelf(lasterr));');
   P = Ppc(:,:);
   
otherwise
   error('Invalid 2nd argument.');
end;

S(1) = [];
if ~isempty(S),
   if (isa(P,'double') | isa(P,'char')) & strcmp(S(1).type,'.'),
      error('Invalid combination in subscripted reference.');
   end;
   eval('P = subsref(P,S);', ...
      ['if isa(P,''pol'') | isa(P,''tsp'') | isa(P,''frac''),', ...
      'error(peel(lasterr)); else error(peelf(lasterr)); end;']);   
end;

%end .. @pol/subsref

