function t = iscoprime(F,tol)
%ISCOPRIME   Test if matrix-den fraction is coprime
%
% The expression  ISCOPRIME(F)  returns 1 if all entries
% of the matrix-den fraction F are coprime, otherwise 0.
%
% An optional argument TOL may specify zeroing tolerance.
% Its default value is the global zeroing tolerance.
%
%        Author:  J.Jezek, 06-Oct-2002
%        Copyright(c) 2002 by Polyx, Ltd.
%        $ Revision $  $ Date 14-Oct-2002 $
%                      $ Date 28-Feb-2003 $

global PGLOBAL;

if nargin==2 & ~isempty(tol),
   if isa(tol,'char'),
      tol = str2num(tol);
   else
      if ~isa(tol,'double'),
         error('Invalid tolerance.');
      end;
   end;
   if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
else
   tol = PGLOBAL.ZEROING;
end;

if isempty(F),
   t = logical(1); return;
end;
if strcmp(F.frac.c,'ncop') & F.frac.tc==tol,
   t = logical(0); return;
end;
if strcmp(F.frac.c,'cop') & F.frac.tc==tol,
   t = logical(1); return;
end;

stol = num2str(tol);
[s1,s2] = size(F); G = pol(zeros(s1,s2));
N = F.frac.num; D = F.frac.den; 
for i = 1:s1,
   for j = 1:s2,
      G(i,j) = gld(N(i,j),D(i,j),stol);
   end;
end;
degG = deg(G,'ent');
t = all(all(degG<=0));
name = inputname(1);
if ~isempty(name),
   if t, props(F,'cop',tol);
   else props(F,'ncop',tol);
   end;
   assignin('caller',inputname(1),F);
end;

%end .. @mdf/iscoprime

