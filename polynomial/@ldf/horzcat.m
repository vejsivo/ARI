function H = horzcat(varargin)
%HORZCAT    Horizontal concatenate left-den fractions
%
% The command
%    H = [F G]   or   H = [F,G]
% returns the horizontal concatenation of left-denominator fractions
% F and G. Both of them must have the same number of rows. Any number
% of matrices may be concatenated within one pair of brackets.
% Horizontal and vertical concatenation may be combined together as
% for standard MATLAB matrices. The command
%    H = HORZCAT(F,G)
% works similarly.
%
% The operation of concatenation is not trivial, it includes
% computation, namely solving a polynomial equation. It uses 
% either the global zeroing tolerance or tolerance TOL, given as
% an optional input argument. To be distinguishable from other
% ones, this argument must have the form of character string, e.g.
%    H = HORZCAT(F,G,'10^-8') .
%
% If the variable symbol is not the same in all matrices, it is set
% to the current value of PGLOBAL.VARIABLE and a warning message is
% issued. However, if there are both 'z' and 'z^-1' among the input
% variables (and no others), the names play a role, no warning
% being issued.
%
% See also LDF/VERTCAT.

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

ni = length(varargin); givtol = [];
if ni>=2,
   lastarg = varargin{ni};
   if isa(lastarg,'char'),
      ni = ni-1; givtol = str2double(lastarg);
      if ~isempty(lastarg),
         if isnan(givtol) | length(givtol)~=1 | ...
               ~isreal(givtol) | givtol<0 | givtol>1,
            error('Invalid last argument.');
         end;
      end;            
   end;
end;
if ~isempty(givtol), tol = givtol;
else tol = PGLOBAL.ZEROING; lastarg = num2str(tol);
end;

allprop = logical(1); pom = 0;
eval('pom = ldf(varargin{1});',...
   'error(''Invalid argument.'');');
Ps = pom.frac.s; Pv = pom.frac.v; inc = 0; Ph = pom.frac.h;
Ptp = pom.frac.tp;
if ~strcmp(pom.frac.p,'prop'), allprop = logical(0);
end;
var = cell(1,ni); var{1} = pom;

if ni>=2,
   for i = 2:ni,
      eval('pom = ldf(varargin{i});',...
         'error(''Invalid argument.'');');
      if pom.frac.s(1)~=Ps(1),
         error('Matrix has not the same number of rows as the other ones.');
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
      var{i} = pom;
   end;
end;

if inc,
   for i = 1:ni,
      v = PGLOBAL.VARIABLE;
      var{i}.frac.v = v;
      var{i}.frac.num.v = v; var{i}.frac.den.v = v;
   end;
end;

if ~isempty(Ph) & Ph==-1,
   warning('Inconsistent sampling periods.');
   for i = 1:ni,
      var{i}.frac.h = [];
      var{i}.frac.num.h = [];
      var{i}.frac.den.h = [];
   end;
end;

if ni==1,
   H = var{1};
   if allprop, props(H,'prop',Ptp);
   end;
   if strcmp(PGLOBAL.COPRIME,'cop'),
      H = coprime(H,tol);
   end;
   if strcmp(PGLOBAL.REDUCE,'red'),
      H = reduce(H,tol);
   else H = smreduce(H);
   end;
   if strcmp(PGLOBAL.DEFRACT,'defr'),
      H = defract(H);
   end;
elseif ni==2,
   X = 0; Y = 0;
   eval('[X Y] = xayb0(var{1}.frac.den,-var{2}.frac.den,tol);', ...
      'error(peel(lasterr));');
   H = ldf(mtimes(X,var{1}.frac.den,tol),...
      [mtimes(X,var{1}.frac.num,tol),mtimes(Y,var{2}.frac.num,tol)]);
   if allprop, props(H,'prop',Ptp);
   end;
   if strcmp(PGLOBAL.COPRIME,'cop'),
      H = coprime(H,tol);
   end;
   if strcmp(PGLOBAL.REDUCE,'red'),
      H = reduce(H,tol);
   else H = smreduce(H);
   end;
   if strcmp(PGLOBAL.DEFRACT,'defr'),
      H = defract(H);
   end;
else
   nihalf = floor(ni/2);
   low = horzcat(var{1:nihalf},lastarg);
   high = horzcat(var{nihalf+1:ni},lastarg);
   H = horzcat(low,high,lastarg);
end;

%end .. @ldf/horzcat

