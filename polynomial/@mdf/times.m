function R = times(P,Q,arg3,arg4);
%TIMES(.*)   Element-wise multiply matrix-den fractions
%           R = P.*Q   or   R = TIMES(P,Q)
%
% For matrix-den fractions P,Q the command returns matrix-den fraction R
% whose (i,j)-th entry is   P(i,j)*Q(i,j)
%
% Matrices P,Q must have the same dimensions, unless one of them is
% scalar. In such a case, this scalar is multiplied with every entry
% of the second matrix.
%
% The variable symbols of P,Q should be the same. When not, a warning
% is issued and the symbols are changed to the standard one. However, 
% if one symbol is 'z' and the second 'z^-1' then the symbols play
% a role. The resulting symbol is taken form P, no warning being issued.
%
% An optional input argument TOL can specify the zeroing tolerance
% to be used instead of the standard one.
%
% The result is standardly of class 'mdf'. However, the class of the
% result may be specified by an optional input argument CLASS, with
% values 'rdf', 'ldf', 'mdf' or 'sdf'.
%
% See also MDF/MTIMES.

%       Author:  J. Jezek  26-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 26-Apr-2000 $
%                     $ Date 06-Nov-2000 $
%                     $ Date 24-Jan-2002 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

tol = PGLOBAL.ZEROING; type = 'mdf';
na = nargin;
if na<2,
   error('Not enough input arguments.');
elseif na>2 & ~isempty(arg3),
   if isa(arg3,'double'), tol = arg3;
   elseif isa(arg3,'char'), type = arg3;
   else error('Invalid 3rd argument.');
   end;
end;
if na>3 & ~isempty(arg4),
   if isa(arg4,'double'), tol = arg4;
   elseif isa(arg4,'char'), type = arg4;
   else error('Invalid 4th argument.');
   end;
end;

cop = PGLOBAL.COPRIME; PGLOBAL.COPRIME = 'cop';
eval('P = mdf(P);','error(''Invalid 1st argument.'');');
eval('Q = mdf(Q);','error(''Invalid 2nd argument.'');');

[tv,Rv,P,Q] = testvf(P,Q);
if ~tv, warning('Inconsistent variables.');
end;
[th,Rh,P,Q] = testhf(P,Q,Rv);
if ~th, warning('Inconsistent sampling periods.');
end;
[td,P,Q] = testdf(P,Q);
if ~td, error('Matrices not of the same dimensions.');
end;

Rn = 0; Rd = 0;
eval('Rn = times(P.frac.num,Q.frac.num,tol); Rd = times(P.frac.den,Q.frac.den,tol);', ...
   'error(peel(lasterr));');
R = mdf(Rn,Rd);
PGLOBAL.COPRIME = cop;

if strcmp(P.frac.p,'prop') & strcmp(Q.frac.p,'prop') & ...
      P.frac.tp==Q.frac.tp,
   props(R,'prop',P.frac.tp);
end;

if strcmp(type,'sdf'), R = sdf(R,tol);
elseif strcmp(type,'rdf'), R = rdf(R,tol);
elseif strcmp(type,'ldf'), R = ldf(R,tol);
elseif ~strcmp(type,'mdf'),
   error('Invalid command option.');
end;

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
  
%end .. @mdf/times
