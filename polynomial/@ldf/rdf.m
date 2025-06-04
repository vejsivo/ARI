function G = rdf(F,tol)
%RDF   Convert left-den fraction to right-den fraction
%
% Given a left-denominator fraction  F = D^-1 * N,  the command
%     G = RF(F)
% computes a right-denominator fraction  G = P * Q^-1 such that
%     G = F,   i.e.   P * Q^-1 = D^-1 * N
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also RDF/RDF, RDF/LDF.

%        Author:  J. Jezek  10-Nov-1999
%        Copyright(c) 1999 by Polyx Ltd.
%        $ Revision $  $ Date 21-Apr-2000 $
%                      $ Date 26-Jul-2000 $
%                      $ Date 02-Nov-2000 $
%                      $ Date 30-Sep-2002 $
%                      $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin==1,
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double')
   error('Invalid tolerance.');
end;

Fd = F.frac.den; [sF1,sF2] = size(F);
if ~strcmp(PGLOBAL.COPRIME,'cop') & ~isempty(F) & ...
      (sF1==1 | all(all(Fd==diag(repmat(Fd(1),sF1,1))))),
   G = rdf(F.frac.num,Fd(1)*eye(sF2));
   props(G,F.frac.p,F.frac.tp,F.frac.r);
else
   PP = 0; QQ = 0;
   eval('[PP,QQ] = axby0(Fd,-F.frac.num,tol);', 'error(peel(lasterr));');
   G = rdf(PP,QQ);
   props(G,'cop',F.frac.p,F.frac.tp);
end;

if strcmp(PGLOBAL.REDUCE,'red'), G = reduce(G,tol);
else G = smreduce(G);
end;

%end .. @ldf/rdf
