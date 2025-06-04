function t = iscoprime(F,tol)
%ISCOPRIME   Test if right-den fraction is coprime
%
% The expression  ISCOPRIME(F)  returns 1 if the right-den
% fraction F is right coprime, otherwise 0.
%
% An optional argument TOL may specify zeroing tolerance.
% Its default value is the global zeroing tolerance.
%
%        Author:  J.Jezek, 06-Oct-2002
%        Copyright(c) 2002 by Polyx, Ltd.
%        $ Revision $ $ Date 14-Oct-2002 $
%                     $ Date 28-Feb-2003 $

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
G = grd(F.frac.num,F.frac.den,stol);
degG = deg(G);
t = isempty(degG) | degG==0;
name = inputname(1);
if ~isempty(name),
   if t, props(F,'cop',tol);
   else props(F,'ncop',tol);
   end;
   assignin('caller',inputname(1),F);
end;

%end .. @rdf/iscoprime

