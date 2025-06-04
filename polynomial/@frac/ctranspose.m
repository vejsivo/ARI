function [Fct,varargout] = ctranspose(F)
%CTRANSPOSE   Conjugated transpose fraction
%
% For F, matrix-den fraction or scalar-den fraction,
% the command  FT = F'  or  FT = CTRANSPOSE(F)
% returns the conjugated transpose of F.
%
% If F is a fraction in variable 's', then
%    FT(s) = CONJ(F'(-s)) .
% For variable 'p' the command works similarly.
%
% If F is a fraction in variable 'z', then
%    FT(z) = CONJ(F'(1/z)) .
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
% See also FRAC/TRANSPOSE, FRAC/CONJ, POL/CTRANSPOSE.

%        Author: J. Jezek  30-Mar-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        $ Revision $  $ Date 06-Jul-2001 $
%                      $ Date 25-Jul-2002 $
%                      $ Date 30-Sep-2002 $
%                      $ Date 14-Oct-2002 $

global PGLOBAL;

no = nargout;
Fct = F;

%variables 'd', 'q'
if strcmp(F.v,'d') | strcmp(F.v,'q'),
   if no > 2,
      error('Too many output arguments.');
   end;
   [Fct.num,nn] = ctranspose(F.num);
   [Fct.den,nd] = ctranspose(F.den);
   if nd>=nn, Fct.num = shift(Fct.num,nd-nn,F.v);
   else       Fct.den = shift(Fct.den,nn-nd,F.v);
   end;
   props(Fct,'prop?',[]);
   varargout{1} = 0;
   
%variables 's', 'p', 'z', 'z^-1'
elseif ~isempty(F.v),
   if no > 1,
      error('Too many output arguments.');
   end;
   Fct.num = ctranspose(F.num);
   Fct.den = ctranspose(F.den);
   if strcmp(F.v,'z'), props(Fct,'z^-1','prop?');
   elseif strcmp(F.v,'z^-1'), props(Fct,'z','prop?');
   else props(Fct,F.p,F.tp);
   end;

%variable ''
else
   if no > 2,
      error('Too many output arguments.');
   end;
   Fct.num = ctranspose(F.num);
   Fct.den = ctranspose(F.den);
   props(Fct,'prop',PGLOBAL.ZEROING);
   varargout{1} = 0;
end;

Fct.s = fliplr(F.s);

Fcl = class(F);
if ~strcmp(Fcl,'frac'),
   if strcmp(PGLOBAL.COPRIME,'cop'), Fct = coprime(Fct);
   end;
   if strcmp(PGLOBAL.REDUCE,'red'), Fct = reduce(Fct);
   else Fct = smreduce(Fct);
   end;
   if strcmp(PGLOBAL.DEFRACT,'defr'), Fct = defract(Fct);
   end;
end;

%end .. @frac/ctranspose

