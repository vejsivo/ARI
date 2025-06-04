function R = subsasgn(R,S,Q)
%SUBSASGN    Subscripted assignment for scalar-den fraction
%
% For scalar-denominator fraction R,
%    R.num = Q         sets the numerator of R
%    R.den = Q         sets the denominator
%    R.v   = Q         sets the variable symbol
%    R.h   = Q         sets the sampling period
%    R.u   = Q         sets the user data
%    R.c   = Q         sets the coprime flag
%    R.r   = Q         sets the reduce flag
%    R.p   = Q         sets the proper flag
%    R.tc  = Q         sets the coprime tolerance
%    R.tp  = Q         sets the proper tolerance
%    R(i,j) = Q        sets the (i,j)-th entry of R
%    R(i1:i2,:) = Q    similarly
%
% For the numerator and denominator, the subscripting in assignment
% may be further elaborated, as for polynomials, e.g.  R.num(i,j) = Q
% or  R.den{k} = Q 
%
% If the variable names are assigned in the numerator and in the
% denominator inconsistently, a warning is issued and the default
% variable name is used.

%        Author:  J. Jezek  24-Jan-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        $ Revision $  $ Date 21-Apr-2000 $
%                      $ Date 26-May-2000 $
%                      $ Date 16-Oct-2000 $
%                      $ Date 05-Feb-2001 $
%                      $ Date 08-Jul-2001 $
%                      $ Date 06-Oct-2002  bug $
%                      $ Date 14-Oct-2002 $

global PGLOBAL;

vvx = PGLOBAL.VARIABLE;
Sl = length(S); St1 = S(1).type;
if strcmp(St1,'.') & strcmp(S(1).subs,'frac'),
   S = S(2:end); Sl = length(S); St1 = S(1).type;
end;
if strcmp(St1,'.'),
   flstr = lower(S(1).subs);
   if strcmp(flstr,'den'),
      if Sl==1,
         eval('R.frac.den = pol(Q);','error(peel(lasterr));');
      else
         eval('R.frac.den = subsasgn(R.frac.den,S(2:end),Q);',...
              'error(peel(lasterr));');
      end;
      eval('R = sdf(R.frac.num,R.frac.den);','error(peel(lasterr));');
      return;  
   end;
   
   switch flstr(1);
   case 'n', if Sl==1,
                eval('R.frac.num = pol(Q);','error(peel(lasterr));');
             else
                eval('R.frac.num = subsasgn(R.frac.num,S(2:end),Q);',...
                     'error(peel(lasterr));');
             end;
             eval('R = sdf(R.frac.num,R.frac.den);','error(peel(lasterr));');
             return;
             
   case 'v', eval('R.frac.v = Q;','error(peel(lasterr));');
   case 'h',
      if isempty(Q),
         I = strmatch(R.frac.v,{'z^-1';'z';'d';'q'},'exact');
         if ~isempty(I),
            R.frac.num.h = []; R.frac.den.h = []; R.frac.h = [];
         else return;
         end;
      elseif ~isa(Q,'double') | length(Q)~=1 | ...
            ~isreal(Q) | Q<0 | isinf(Q),
         error('Invalid sampling period.');
      elseif strcmp(R.frac.v,'s') | strcmp(R.frac.v,'p'),
         if Q>0,
            error('Invalid object for nonzero sampling period.');
         end;
         R.frac.num.h = 0; R.frac.den.h = 0; R.frac.h = 0;
      elseif isempty(R.frac.v),
         R.frac.num.h = []; R.frac.den.h = []; R.frac.h = [];
      else
         if Q==0,
            switch R.frac.v;
            case 'z', RR = R;
            case 'q', RR = R; RR.frac.v = 'z';
            case 'z^-1', RR = reverse(R);
            case 'd', RR = R; RR.frac.v = 'z^-1'; RR = reverse(RR);
            end;
            num = shift(RR.frac.num,-1);
            if isa(num,'tsp'),
               error('Invalid object for zero sampling period.');
            end;
            zz = z(num.h);
            num = num*(zz-1); den = RR.frac.den;
            [r,c] = size(RR);
            quot = pol(ones(r,c)); rem = quot;
            for i = 1:r,
               for j = 1:c,                
                  [quot(i,j),rem(i,j)] = rdiv(num(i,j),den);
               end;
            end;
            if deg(rem)>=0 | deg(quot)>0,
               error('Invalid object for zero sampling period.');
            end;
            R.frac.nm.h = []; R.frac.den.h = []; R.frac.h = 0;
         else
            R.frac.num.h = Q; R.frac.den.h = Q; R.frac.h = Q;
         end;
      end;
      
   case 'c', eval('R.frac.c = Q;','error(peel(lasterr));');
   case 'r', eval('R.frac.r = Q;','error(peel(lasterr));');
   case 'p', eval('R.frac.p = Q;','error(peel(lasterr));');
   case 't',
      switch flstr(2),
      case 'c', eval('R.frac.tc = Q;','error(peel(lasterr));');
      case 'p', eval('R.frac.tp = Q;','error(peel(lasterr));');
      otherwise error('Invalid field in subscripted assignment.');
      end;
      
   case 'u',
      S(1) = [];
      if isempty(S), R.frac.u = Q;
      else eval('R.frac.u = subsasgn(R.frac.u,S,Q);', ...
            ['if isa(R.frac.u,''pol'') | isa(R.frac.u,''tsp'') | isa(R.frac.u,''frac''),', ...
               'error(peel(lasterr)); else error(peelf(lasterr)); end;']);
      end;
      return;
      
   otherwise error('Invalid field in subscripted assignment.');
   end;
   if Sl>1,
      error('Invalid subscripted assignment.');
   end;
   
elseif strcmp(St1,'()'),
   Ssubs = S(1).subs; lS = length(Ssubs);
   if lS==1 | lS==2,
      eval('Q = sdf(Q);','error(peel(lasterr));');
      [X Y] = axby0(R.frac.den,-Q.frac.den);
      R.frac.num = R.frac.num*X; R.frac.den = R.frac.den*X;
      Q.frac.num = Q.frac.num*Y;
      eval('R.frac.num = subsasgn(R.frac.num,S,Q.frac.num);','error(peel(lasterr));');
      eval('R = sdf(R.frac.num,R.frac.den);','error(peel(lasterr));');    
   else
      error('More than two subscripts.');
   end;
else
   error('Invalid subscripted assignment.');
end;

%end .. @sdf/subsasgn
