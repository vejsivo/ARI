function R = subsasgn(R,S,Q)
%SUBSASGN    Subscripted assignment for fraction
%
% For scalar fraction R,
%    R.num = Q         sets the numerator of R
%    R.den = Q         sets the denominator
%    R.s   = Q         sets the size
%    R.v   = Q         sets the variable symbol
%    R.h   = Q         sets the sampling period
%    R.u   = Q         sets the user data
%    R.c   = Q         sets the coprime flag
%    R.r   = Q         sets the reduce flag
%    R.p   = Q         sets the proper flag
%    R.tc  = Q         sets the coprime tolerance
%    R.tp  = Q         sets the proper tolerance
% For the numerator and denominator, the subscripting in assignment
% may be further elaborated, as for polynomials, e.g.  R.num(i,j) = Q
% or  R.den{k} = Q 
%
% If the variable names are assigned in the numerator and in the
% denominator inconsistently, a warning is issued and the default
% variable name is used.

%        Author:  J. Jezek  24-Jan-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        $ Revision $  $ Date 14-Oct-2002 $

global PGLOBAL;

vvx = PGLOBAL.VARIABLE;
Sl = length(S); St1 = S(1).type;
if strcmp(St1,'.'),
   flstr = lower(S(1).subs);
   if strcmp(flstr,'den'),
      if Sl==1,
         eval('R.den = pol(Q);','error(peel(lasterr));');
      else
         eval('R.den = subsasgn(R.den,S(2:end),Q);',...
            'error(peel(lasterr));');
      end;
   elseif strcmp(flstr,'tp'),
      if Sl>1,
          error('Invalid subscripted assignment.');
      end;
      if ~strcmp(R.p,'prop?'),
         if ~isa(Q,'double') | ~length(Q)~=1 | ~isreal(Q) | ...
                 Q<0 | Q>1,
             error('Invalid tolerance.');
         end;
         R.tp = Q;
      end;
   elseif strcmp(flstr,'tc'),
      if Sl>1,
         error('Invalid subscripted assignment.');
      end;
      if ~strcmp(R.c,'cop?'),
         if ~isa(Q,'double') | length(Q)~=1 | ~isreal(Q) | ...
                 Q<0 | Q>1,
             error('Invalid tolerance.');
         end;
         R.tc = Q;
      end;       
   else
      switch flstr(1),
      case 'n',
         if Sl==1,
            eval('R.num = pol(Q);','error(peel(lasterr));');
         else
            eval('R.num = subsasgn(R.num,S(2:end),Q);',...
               'error(peel(lasterr));');
         end;
      case 's',
         if Sl>1,
            error('Invalid subscripted asignment.');
         end;
         R.s = Q;
      case 'v',
         if Sl>1,
            error('Invalid subscripted assignment.');
         end;
         if isa(Q,'pol'),
            [vs1,vs2,vd] = size(Q);
            if all([vs1,vs2,vd]==1)&(~any(Q.c(:,:)-[0,1])),
               Q = Q.v;
            else
               error('Invalid variable symbol.');
            end;
         end;
         if isstr(Q),
            if strcmp(Q,'zi'), Q = 'z^-1';
            end;
            I = strmatch(Q,{'s';'p';'z^-1';'d';'z';'q'},'exact');
            if ~isempty(I) | (isempty(Q) & ...
                  ((isempty(R.num)|R.num.d<=0) & ...
                  (isempty(R.den)|R.den.d==0)) ),
               Oldv = R.v; R.v = Q; R.num.v = Q; R.den.v = Q;
               J = strmatch(Q,{'z^-1';'d';'z';'q'},'exact');
               if ~isempty(J) & (strcmp(Oldv,'s') | strcmp(Oldv,'p')),
                  R.h = 1;
               end;   
            else,
               error('Invalid variable symbol.');
            end;
            if isempty(R.num.v) & isempty(R.den.v),
               R.v = '';
            end;
         else, 
            error('Invalid variable symbol.'); 
         end;
         R.c = 'cop?'; R.r = 'red?'; R.p = 'prop?';
      case 'h',
         if Sl>1,
            error('Invalid subscripted assignment.');
         end;
         if isempty(Q), Q = [];
         end;
         R.h = Q;
      case 'c',
         if Sl>1,
            error('Invalid subscripted assignment.');
         end;
         if isa(Q,'char'),
            if strcmp(Q,'cop?') | strcmp(Q,'coprime?'),
               R.c = 'cop?'; R.tc = [];
            else
               error('Invalid flag; only ''cop?'' allowed.');
            end;
         else
            error('Invalid flag; only ''cop?'' allowed.');
         end;
      case 'r',
         if Sl>1,
            error('Invalid subscripted assignment.');
         end;
         if isa(Q,'char'),
            if strcmp(Q,'red?') | strcmp(Q,'reduce?'),
               R.r = 'red?';
            else
               error('Invalid flag; only ''red?'' allowed.');
            end;
         else
            error('Invalid flag; only ''red?'' allowed.');
         end;
      case 'p',
         if Sl>1,
            error('Invalid subscripted assignment.');
         end;
         if isa(Q,'char'),
            if strcmp(Q,'prop?') | strcmp(Q,'proper?'),
               R.p = 'prop?'; R.tp = [];
            else
               error('Invalid flag; only ''prop?'' allowed.');
            end;
         else
            error('Invalid flag; only ''prop?'' allowed.');
         end;
      case 'u',
          if Sl==1,
              R.u = Q;
          else
              eval('R.u = subsasgn(R.u,S(@:end),Q);', ...
                  'error(peel(lasterr));');
          end;
      otherwise error('Invalid field in subscripted assignment.');
      end;
   end;
else
   error('Invalid subscripted assignment.');
end;

%end .. @frac/subsasgn
