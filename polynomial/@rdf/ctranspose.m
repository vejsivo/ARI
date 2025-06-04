function [Fct,varargout] = ctranspose(F)
%CTRANSPOSE   Conjugated transpose right-den fraction
%
% The command  FT = F'  or  FT = CTRANSPOSE(F)
% returns the conjugated transpose of the right-den fraction F.
% The result is left-den fraction.
%
% If F is a fraction in variable 's', then
%    FT(s) = CONJ(F'(-s)) .
% For variable 'p' the command works similarly.
%
% If F is a fraction in variable 'z', then
%    FT(z) = CONJ(F'(z^-1)) .
% The result is a fraction in 'z^-1'.
% For variable 'z^-1' the command works similarly,
% the result is a fraction in 'z'.
%
% If F is a fraction in variable 'q', then
%    FT(q) = CONJ(F'(q^-1)) .
% The result is again a fraction in 'q'. The command
%    [FT,N] = CTRANSPOSE(F)
% returns both FT and N = 0. This is for compatibility
% with POL/CTRANSPOSE.
% For variable 'd' the command works similarly.
%
% See also RDF/TRANSPOSE, FRAC/CONJ, POL/CTRANSPOSE.

%        Author: J. Jezek  30-Mar-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        $ Revision $  $ Date 26-Apr-2000 $
%                      $ Date 06-Jul-2001 $
%                      $ Date 30-Sep-2002 $
%                      $ Date 14-Oct-2002 $

global PGLOBAL;

no = nargout;
N = F.frac.num; D = F.frac.den;

%variables 'd', 'q'
if strcmp(F.frac.v,'d') | strcmp(F.frac.v,'q'),
   if no > 2,
      error('Too many output arguments.');
   end;
   [N,nn] = ctranspose(N);
   [D,nd] = ctranspose(D);
   if nd>=nn, N = shift(N,nd-nn,F.frac.v);
   else       D = shift(D,nn-nd,F.frac.v);
   end;
   Fct = ldf(D,N);
   varargout{1} = 0;
   
%variables 's', 'p', 'z', 'z^-1'
elseif ~isempty(F.frac.v),
   if no > 1,
      error('Too many output arguments.');
   end;
   N = ctranspose(N);
   D = ctranspose(D);
   Fct = ldf(D,N);
   if strcmp(F.frac.v,'s') | strcmp(F.frac.v,'p'),
      props(Fct,F.frac.p,F.frac.tp);
   end;
   
%variable ''
else
   if no > 2,
      error('Too many output arguments.');
   end;
   N = ctranspose(N);
   D = ctranspose(D);
   Fct = ldf(D,N);
   props(Fct,'prop',PGLOBAL.ZEROING);
   varargout{1} = 0;
end;

if strcmp(PGLOBAL.COPRIME,'cop'), Fct = coprime(Fct);
end;
if strcmp(PGLOBAL.REDUCE,'red'), Fct = reduce(Fct);
else Fct = smreduce(Fct);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), Fct = defract(Fct);
end;

%end .. @rdf/ctranspose

