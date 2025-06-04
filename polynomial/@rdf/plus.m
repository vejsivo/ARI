function H = plus(F,G,arg3,arg4);
%PLUS (+)   Add right-den fractions
%            C = A + B
%
% C = A + B or C = PLUS(A,B) adds the right-denominator fractions A and B 
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
% be used. Standardly, the result is right fraction. However, an optional
% 'ldf' or 'rdf' input argument may specify the class of the result.
%
% See also RDF/MINUS, RDF/UPLUS.

%       Author:  J. Jezek  18-Nov-1999
%       Copyright(c) 1999 by Polyx Ltd.
%       $ Revision $  $ Date 25-Apr-2000 $
%                     $ Date 30-May-2000 $
%                     $ Date 02-Aug-2000 $
%                     $ Date 01-Nov-2000 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

tol = PGLOBAL.ZEROING; hand = 'rdf';
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

eval('F = rdf(F);','error(peel(lasterr));');

if strcmp(hand,'rdf'),
   if isa(G,'ldf'),
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
      GF = plus(mtimes(Gd,Fn,tol),mtimes(Gn,Fd,tol),tol);
      [X Y] = axby0(Gd,-GF,tol);
      Hd = mtimes(Fd,Y,tol);
      H = rdf(X,Hd);
      if strcmp(F.frac.p,'prop') & strcmp(G.p,'prop') & ...
            F.frac.tp==G.tp,
         props(H,'prop',F.frac.tp);
      end;
   else
      eval('G = rdf(G);','error(peel(lasterr));');
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
      Gn = G.frac.num; Gd = G.frac.den; 
      X = 0; Y = 0;
      eval('[X,Y] = axby0(Fd,-Gd,tol);', 'error(peel(lasterr));');
      FnX = mtimes(Fn,X,tol); GnY = mtimes(Gn,Y,tol);
      Hn = plus(FnX,GnY,tol); Hd = mtimes(Fd,X,tol);
      H = rdf(Hn,Hd);
      if strcmp(F.frac.p,'prop') & strcmp(G.frac.p,'prop') & ...
            F.frac.tp==G.frac.tp,
         props(H,'prop',F.frac.tp);
      end;
   end;
elseif strcmp(hand,'ldf'),
   eval('G = ldf(G);','error(peel(lasterr));');
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
   GF = plus(mtimes(Gd,Fn,tol),mtimes(Gn,Fd,tol),tol);
   [X Y] = xayb0(Fd,-GF,tol);
   Hd = mtimes(Y,Gd,tol);
   H = ldf(Hd,X);
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

%end .. @rdf/plus   
