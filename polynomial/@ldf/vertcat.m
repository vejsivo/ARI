function H = vertcat(varargin)
%VERTCAT    Vertical concatenate left-den fractions
%
% The command
%    H = [F;G]
% returns the vertical concatenation of left-denominator fractions
% F and G. Both of them must have the same number of columns. Any
% number of matrices may be concatenated within one pair of brackets.
% Vertical and horizontal concatenation may be combined together as
% for standard MATLAB matrices. The command
%    H = VERTCAT(F,G)
% works similarly.
%
% If the variable symbol is not the same in all matrices, it is set
% to the current value of PGLOBAL.VARIABLE and a warning message is
% issued. However, if there are both 'z' and 'z^-1' among the input
% variables (and no others), the names play a role, no warning
% being issued.
%
% See also LDF/HORZCAT.

%         Author:  J. Jezek  30-Nov-1999
%         Copyright(c) 1999 by Polyx, Ltd.
%         $ Revision $  $ Date 21-Apr-2000 $
%                       $ Date 30-May-2000 $
%                       $ Date 03-Aug-2000 $
%                       $ Date 02-Nov-2000 $
%                       $ Date 29-Jan-2002 $
%                       $ Date 30-Sep-2002 $
%                       $ Date 14-Oct-2002 $
%                       $ Date 28-Feb-2003 $

global PGLOBAL;

ni = length(varargin); tol = [];
if ni>=2,
   lastarg = varargin{ni};
   if isa(lastarg,'char'),
      ni = ni-1; tol = str2double(lastarg);
      if ~isempty(lastarg),
         if isnan(tol) | length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
            error('Invalid last argument.');
         end;
      end;      
   end;
end;

allprop = logical(1); pom = 0;
eval('pom = ldf(varargin{1});',...
   'error(''Invalid argument.'');');
Ps = pom.frac.s; Pv = pom.frac.v; inc = 0; Ph = pom.frac.h;
Ptp = pom.frac.tp;
if ~strcmp(pom.frac.p,'prop'), allprop = logical(0);
end;
num = cell(1,ni); den = cell(1,ni);
num{1} = pom.frac.num; den{1} = pom.frac.den;

if ni>=2,
   for i = 2:ni,
      eval('pom = ldf(varargin{i});',...
         'error([''Invalid argument.'');');
      if pom.frac.s(2)~=Ps(2),
         error('Matrix has not the same number of columns as the other ones.');
      end;
      if allprop,
         if ~strcmp(pom.frac.p,'prop') | pom.frac.tp~=Ptp,
            allprop = logical(0);
         end;
      end;
      if ~isempty(pom.frac.v),
         if isempty(Pv),
            Pv = pom.frac.v;
         elseif (strcmp(Pv,'z') & strcmp(pom.frac.v,'z^-1')) | ...
                (strcmp(Pv,'z^-1') & strcmp(pom.frac.v,'z')),
            pom = reverse(pom);
         elseif ~strcmp(Pv,pom.frac.v),
            if ~inc,
               warning('Inconsistent variables.');
               inc = 1;
            end;
         end;
      end;
      if isempty(Ph),
         Ph = pom.frac.h;
      elseif Ph~=-1,
         if ~isempty(pom.frac.h) & pom.frac.h~=Ph,
            Ph = -1;
         end;
      end;   
      num{i} = pom.frac.num; den{i} = pom.frac.den;
   end;
end;

if inc,
   for i = 1:ni,
      num{i}.v = PGLOBAL.VARIABLE;
      den{i}.v = PGLOBAL.VARIABLE;
   end;
end;

if ~isempty(Ph) & Ph==-1,
   warning('Inconsistent sampling periods.');
   for i = 1:ni,
      num{i}.h = [];
      den{i}.h = [];
   end;
end;

H = ldf(blkdiag(den{1:ni}),vertcat(num{1:ni}));

if allprop, props(H,'prop',Ptp);
end;
if strcmp(PGLOBAL.COPRIME,'cop'),
   H = coprime(H,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   H = reduce(H,tol);
else
   H = smreduce(H);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'),
   H = defract(H);
end;

%end .. @ldf/vertcat
