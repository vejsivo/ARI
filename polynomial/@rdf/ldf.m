function G = ldf(F,tol)
%LDF   Convert right-den fraction to left-den fraction
%
% Given a right-denominator fraction  F = N * D^-1, the command
%     G = LDF(F)
% computes a left_denominator fraction  G = Q^-1 * P  such that
%     G == F     i.e     Q^-1 * P == N * D^-1
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also LDF/LDF, LDF/RDF.

%       Author:  J. Jezek  10-Nov-1999
%       Copyright(c) 1999 by Polyx Ltd.
%       $ Revision $  $Date 21-Apr-2000 $
%                     $Date 26-Jul-2000 $
%                     $Date 01-Nov-2000 $
%                     $Date 30-Sep-2002 $
%                     $Date 14-Oct-2002 $

global PGLOBAL;

if nargin==1,
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double'),
   error('Invalid tolerance.');
end;

Fd = F.frac.den; [sF1,sF2] = size(F);
if ~strcmp(PGLOBAL.COPRIME,'cop') & ~isempty(F) & ...
      (sF2==1 | all(all(Fd==diag(repmat(Fd(1),sF2,1))))),
   G = ldf(Fd(1)*eye(sF1),F.frac.num);
   props(G,F.frac.p,F.frac.tp,F.frac.r);
else
   PP = 0; QQ = 0;
   eval('[PP,QQ] = xayb0(Fd,-F.frac.num,tol);', 'error(peel(lasterr));');
   G = ldf(QQ,PP);
   props(G,'cop',tol,F.frac.p,F.frac.tp);
end;

if strcmp(PGLOBAL.REDUCE,'red'), G = reduce(G,tol);
else G = smreduce(G);
end;

%end .. @rdf/ldf
