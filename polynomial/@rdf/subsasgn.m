function F = subsasgn(F,S,R)
%SUBSASGN    Subscripted assignment for right-den fraction
%
% For right-denominator fraction F,
%    F.num = R   sets the numerator of F
%    F.den = R   sets the denominator
%    F.v   = R   sets the variable symbol
%    F.h   = R   sets the sampling period
%    F.u   = R   sets the user data
%    F.c   = R   sets the coprime flag
%    F.r   = R   sets the reduced flag
%    F.p   = R   sets the proper flag
%    F.tc  = R   sets the coprime tolerance
%    F.tp  = R   sets the proper tolerance
%    F(i,j) = R      sets the (i,j)-th entry of F
%    F(i1:i2,:) = R  similarly

%
% For the numerator and denominator, the subscripting in assignment
% may be further elaborated, as for polynomials, e.g.  F.num(i,j) = R
% or  F.den{k} = R 
%
% If the variable names are assigned in the numerator and in the
% denominator inconsistently, a warning is issued and the default
% variable name is used.

%        Author:  J. Jezek  25-Nov-1999
%        Copyright(c) 1999 by Polyx, Ltd.
%        $ Revision $  $ Date 21-Apr-2000 $
%                      $ Date 22-May-2000 $
%                      $ Date 16-Oct-2000 $
%                      $ Date 06-Feb-2001 $
%                      $ Date 08-Jul-2001 $
%                      $ Date 28-Feb-2002 $
%                      $ Date 14-Oct-2002 $

global PGLOBAL;

vvx = PGLOBAL.VARIABLE; erfg =0;
Sl = length(S); St1 = S(1).type;
if strcmp(St1,'.') & strcmp(S(1).subs,'frac'),
   S = S(2:end); Sl = length(S); St1 = S(1).type;
end;
if strcmp(St1,'.'),
   flstr = lower(S(1).subs);
   if strcmp(flstr,'den'),
      if Sl==1,
         eval('F.frac.den = pol(R);','error(peel(lasterr));');
      else
         eval('F.frac.den = subsasgn(F.frac.den,S(2:end),R);',...
              'error(peel(lasterr));');
      end;
      eval('F = rdf(F.frac.num,F.frac.den);','error(peel(lasterr));');
      return;
   end;
   
   switch flstr(1);
   case 'n',
      if Sl==1,
         eval('F.frac.num = pol(R);','error(peel(lasterr));');
      else
         eval('F.frac.num = subsasgn(F.frac.num,S(2:end),R);',...
              'error(peel(lasterr));');
      end;
      Fu = F.frac.u;  
      eval('F = rdf(F.frac.num,F.frac.den);','error(peel(lasterr));'); 
      F.frac.u = Fu;
      return;
             
   case 'v',
      eval('F.frac.v = R;','error(peel(lasterr));');
      
   case 'h',
      if isempty(R),
         I = strmatch(F.frac.v,{'z^-1';'z';'d';'q'},'exact');
         if ~isempty(I),
            F.frac.num.h = []; F.frac.den.h = []; F.frac.h = [];
         else return;
         end;
      elseif ~isa(R,'double') | length(R)~=1 | ...
            ~isreal(R) | R<0 | isinf(R),
         error('Invalid sampling period.');
      elseif strcmp(F.frac.v,'s') | strcmp(F.frac.v,'p'),
         if R>0,
            error('Invalid object for nonzero sampling period.');
         end;
         F.frac.num.h = 0; F.frac.den.h = 0; F.frac.h = 0;
      elseif isempty(F.frac.v),
         F.frac.num.h = []; F.frac.den.h = []; F.frac.h = [];
      else
         if R==0,
            switch F.frac.v;
            case 'z', FF = F;
            case 'q', FF = F; FF.frac.v = 'z';
            case 'z^-1', FF = reverse(F);
            case 'd', FF = F; FF.frac.v = 'z^-1'; FF = reverse(FF);
            end;
            num = shift(FF.frac.num,-1);
            if isa(num,'tsp'),
               error('Invalid object for zero sampling period.');
            end;
            zz = z(num.h);
            num = num*(zz-1); den = FF.frac.den;
            [quot,rem] = rdiv(num,den);
            if deg(rem)>=0 | deg(quot)>0,
               error('Invalid object for zero sampling period.');
            end;
            F.frac.num.h = []; F.frac.den.h = []; F.frac.h = 0;
         else
            F.frac.num.h = R; F.frac.den.h = R; F.frac.h = R;
         end;
      end;
      
   case 'u',
      S(1) = [];
      if isempty(S), F.frac.u = R;
      else eval('F.frac.u = subsasgn(F.frac.u,S,R);', ...
            ['if isa(F.frac.u,''pol'') | isa(F.frac.u,''tsp'') | isa(F.frac.u,''frac''),', ...
               'error(peel(lasterr)); else error(peelf(lasterr)); end;']);
      end;
      return;
      
   case 'c', eval('F.frac.c = R;','error(peel(lasterr));');
   case 'r', eval('F.frac.r = R;','error(peel(lasterr));');
   case 'p', eval('F.frac.p = R;','error(peel(lasterr));'); 
   case 't',
      switch flstr(2),
      case 'c', eval('F.frac.tc = R;','error(peel(lasterr));');
      case 'p', eval('F.frac.tp = R;','error(peel(lasterr));');
      otherwise error('Invalid field in subscripted assignment.');
      end;      
   otherwise error('Invalid field in subscripted assignment.');
   end;
   
elseif strcmp(St1,'()'),
   F = sdf(F);
   eval('F = subsasgn(F,S,R);','error(peel(lasterr));');
   F = rdf(F);
else
   error('Invalid subscripted assignment.');
end;

%end .. @rdf/subsasgn
