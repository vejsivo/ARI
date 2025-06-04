function H = plus(F,G,arg3,arg4);
%PLUS (+)   Add left-den fractions
%            C = A + B
%
% C = A + B or C = PLUS(A,B) adds the left-den fractions A and B 
% with zeroing according to the global zeroing tolerance. The dimensions
% of matrices A and B must be the same unless one of them is scalar;
% in such a case the scalar is added to every entry of the other matrix.
%
% The variable symbols of F and G should be the same. When not,
% a warning is issued and both symbols are changed to the standard one.
% However, if one symbol is 'z' and the other 'z^-1' then the symbols
% play a role, no warning, the resulting symbol is taken from F.
%
% An optional numerical input argument may specify the tolerance to
% be used. Standardly, the result is left-den fraction. However, an optional
% 'rdf' or 'ldf' input argument may specify the class of the result.
%
% See also LDF/MINUS, LDF/UPLUS.

%       Author:  J. Jezek  18-Nov-1999
%       Copyright(c) 1999 by Polyx Ltd.
%       $ Revision $  $ Date 21-Apr-2000 $
%                     $ Date 30-May-2000 $
%                     $ Date 02-Nov-2000 $
%                     $ Date 24-Jan-2002 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

tol = PGLOBAL.ZEROING; hand = 'ldf';
na = nargin;
if na<2,
   error('Not enough input arguments.');
elseif na>2 & ~isempty(arg3),
   if isa(arg3,'double'), tol = arg3;
   elseif isa(arg3,'char'), hand = arg3;
   else error('Invalid 3rd argument.');
   end;
end;
if na>3 & ~isempty(arg4),
   if isa(arg4,'double'), tol = arg4;
   elseif isa(arg4,'char'), hand = arg4;
   else error('Invalid 4th argument.');
   end;
end;

eval('F = ldf(F);','error(''Invalid 1st argument.'');');

if strcmp(hand,'ldf'),
   if isa(G,'rdf'),
      [tv,Hv,F,G] = testvf(F,G);
      if ~tv, warning('Inconsistent variables.');
      end;
      [th,Hh,F,G] = testhf(F,G,Hv);
      if ~th, warning('Inconsistent sampling periods.');
      end;
      [td,F,G] = testdf(F,G);
      if ~td, error('Matrices not of the same dimensions.');
      end;
      Fn = F.frac.num; Fd = F.frac.den;
      Gn = G.num; Gd = G.den; X = 0; Y = 0;
      FG = plus(mtimes(Fn,Gd,tol),mtimes(Fd,Gn,tol),tol);
      eval('[X Y] = xayb0(Gd,-FG,tol);', 'error(peel(lasterr));');
      Hd = mtimes(Y,Fd,tol);
      H = ldf(Hd,X);
      if strcmp(F.frac.p,'prop') & strcmp(G.p,'prop') & ...
            F.frac.tp==G.tp,
         props(H,'prop',F.frac.tp);
      end;
   else
      eval('G = ldf(G);','error(''Invalid 2nd argument.'');');
      [tv,Hv,F,G] = testvf(F,G);
      if ~tv, warning('Inconsistent variables.');
      end;
      [th,Hh,F,G] = testhf(F,G,Hv);
      if ~th, warning('Inconsistent sampling periods.');
      end;
      [td,F,G] = testdf(F,G);
      if  ~td, error('Matrices not of the same dimensions.');
      end;
      Fn = F.frac.num; Fd = F.frac.den;
      Gn = G.frac.num; Gd = G.frac.den; X = 0; Y = 0;     
      eval('[X,Y] = xayb0(Fd,-Gd,tol);', 'error(peel(lasterr));');
      XFn = mtimes(X,Fn,tol); YGn = mtimes(Y,Gn,tol);
      Hn = plus(XFn,YGn,tol); Hd = mtimes(X,Fd,tol);
      H = ldf(Hd,Hn);
      if strcmp(F.frac.p,'prop') & strcmp(G.frac.p,'prop') & ...
            F.frac.tp==G.frac.tp,
         props(H,'prop',F.frac.tp);
      end;
   end;
elseif strcmp(hand,'rdf'),
   eval('G = rdf(G);','error(''Invalid 2nd argument.'');');
   [tv,Hv,F,G] = testvf(F,G);
   if ~tv, warning('Inconsistent variables.');
   end;
   [th,Hh,F,G] = testhf(F,G,Hv);
   if ~th, warning('Inconsistent sampling periods.');
   end;
   [td,F,G] = testdf(F,G);
   if ~td, error('Matrices not of the same dimensions.');
   end;
   Fn = F.frac.num; Fd = F.frac.den;
   Gn = G.num; Gd = G.den;   
   FG = plus(mtimes(Fn,Gd,tol),mtimes(Fd,Gn,tol),tol);
   [X Y] = axby0(Fd,-FG,tol);
   Hd = mtimes(Gd,Y,tol);
   H = rdf(X,Hd);
   if strcmp(F.frac.p,'prop') & strcmp(G.p,'prop') & ...
         F.frac.tp==G.tp,
      props(H,'prop',F.frac.tp);
   end;
else
   error('Invalid command option.');
end;

if strcmp(PGLOBAL.COPRIME,'cop'), H = coprime(H,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'), H = reduce(H,tol);
else H = smreduce(H);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), H = defract(H);
end;

%end .. @ldf/plus   
   