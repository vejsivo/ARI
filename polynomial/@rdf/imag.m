function [FI,FR] = imag(F,tol)
%IMAG   Imaginary part of right-den fraction
%
% For a right-denominator fraction F, the command
% FI = IMAG(F)  returns the imaginary part of F. The command
% [FI,FR] = IMAG(F) returns also the real part.
%
% The computation is not trivial, it includes solving of a polynomial
% equation. It uses either the global zeroing tolerance or the
% tolerance TOL, given as an optional input argument.
%
% See also RDF/REAL, RDF/ISREAL.

%           Author:  J. Jezek 28-Nov-1999
%           Copyright(c( 1999 by Polyx, Ltd.
%           $ Revision $  $ Date 25-Apr-2000 $
%                         $ Date 30-May-2000 $
%                         $ Date 01-Nov-2000 $
%                         $ Date 21-Sep-2001 $
%                         $ Date 30-Sep-2002 $
%                         $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin==2,
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
else
   tol = PGLOBAL.ZEROING;
end;

if isreal(F.frac.num) & isreal(F.frac.den),
   FI = rdf(zeros(size(F)));
else
   A = real(F.frac.num); B = imag(F.frac.num);
   C = real(F.frac.den); D = imag(F.frac.den);
   E = 0; G = 0;
   eval('[E,G] = axby0(D,C,tol);','error(peel(lasterr));');
   CDEG = minus(mtimes(C,E,tol),mtimes(D,G,tol),tol);
   FI = rdf(plus(mtimes(B,E,tol),mtimes(A,G,tol),tol),CDEG);
   if strcmp(F.frac.p,'prop'), props(FI,'prop',F.frac.tp);
   end;
end;

if strcmp(PGLOBAL.COPRIME,'cop'), FI = coprime(FI,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'), FI = reduce(FI,tol);
else FI = smreduce(FI);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'); FI = defract(FI);
end;

if nargout==2,
   if isreal(F.frac.num) & isreal(F.frac.den),
      FR = F;
   else
      FR = rdf(minus(mtimes(A,E,tol),mtimes(B,G,tol),tol),CDEG);
      if strcmp(F.frac.p,'prop'), props(FR,'prop',F.frac.tp);
      end;
   end;
   
   if strcmp(PGLOBAL.COPRIME,'cop'), FR = coprime(FR,tol);
   end;
   if strcmp(PGLOBAL.REDUCE,'red'), FR = reduce(FR,tol);
   else FR = smreduce(FR);   
   end;
   if strcmp(PGLOBAL.DEFRACT,'defr'), FR = defract(FR);
   end;
end;
      
%end .. @rdf/imag
