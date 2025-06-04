function R = horzcat(varargin)
%HORZCAT    Horizontal concatenate matrix-den fractions
%
% The command
%    R = [P Q]   or   R = [P,Q]
% returns the horizontal concatenation of matrix-den fractions
% P and Q. Both of them must have the same number of rows. Any number
% of matrices may be concatenated within one pair of brackets.
% Horizontal and vertical concatenation may be combined together as for
% standard MATLAB matrices. The command
%    R = HORZCAT(P,Q)
% works similarly.
%
% If the variable symbol is not the same in all matrices, it is set
% to the current value of PGLOBAL.VARIABLE and a warning message is
% issued. However, if there are both 'z' and 'z^-1' among the input
% variables (and no others), the names play a role, no warning
% being issued.
%
% See also MDF/VERTCAT.

%         Author:  J. Jezek  12-Jan-2000
%         Copyright(c) 2000 by Polyx, Ltd.
%         $ Revision $  $ Date 26-Apr-2000 $
%                       $ Date 29-May-2000 $
%                       $ Date 03-Aug-2000 $
%                       $ Date 06-Nov-2000 $
%                       $ Data 29-Jan-2002 $
%                       $ Data 30-Sep-2002 $
%                       $ Date 14-Oct-2002 $

global PGLOBAL;

ni = length(varargin); tol = [];
if ni>=2,
   lastarg = varargin{ni};
   if isa(lastarg,'char'),
      ni = ni-1; tol = str2num(lastarg);
      if ~isempty(lastarg),
         if isempty(tol) | length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
            error('Invalid last argument.');
         end;
      end;
   end;
end;

pom = 0;
eval('pom = mdf(varargin{1});',...
   'error(''Invalid argument.'');');
Ps = pom.frac.s; Pv = pom.frac.v; inc = 0; Ph = pom.frac.h;
num = cell(1,ni); den = cell(1,ni);
num{1} = pom.frac.num; den{1} = pom.frac.den;
Rp = pom.frac.p; Rc = pom.frac.c; Rr = pom.frac.r;
Rtp = pom.frac.tp; Rtc = pom.frac.tc;

if ni>=2,
   for i = 2:ni,
      eval('pom = mdf(varargin{i});',...
         'error(''Invalid argument.'');');
      if pom.frac.s(1)~=Ps(1),
         error('Matrix has not the same number of rows as the other ones.');
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
      if ~strcmp(Rp,'prop?')
         if strcmp(pom.frac.p,'prop?') | pom.frac.tp~=Rtp,
            Rp = 'prop?'; Rtp = [];
         elseif strcmp(pom.frac.p,'nprop'),
            Rp = 'nprop';
         end;
      end;
      if ~strcmp(Rc,'cop?')
         if strcmp(pom.frac.c,'cop?') | pom.frac.tc~=Rtc,
            Rc = 'cop?'; Rtc = [];
         elseif strcmp(pom.frac.c,'ncop'),
            Rc = 'ncop';
         end;
      end;
      if ~strcmp(Rr,'red?') & ~strcmp(pom.frac.r,'red'),
         Rr = pom.frac.r;
      end;
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

R = mdf(horzcat(num{1:ni}),horzcat(den{1:ni}));

props(R,Rp,Rtp,Rc,Rtc,Rr);
if strcmp(PGLOBAL.COPRIME,'cop'),
   R = coprime(R,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   R = reduce(R);
else
   R = smreduce(R);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'),
   R = defract(R);
end;

%end .. @mdf/horzcat
