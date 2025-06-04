function R = times(PP,QQ,arg3,arg4);
%TIMES(.*)   Element-wise multiply scalar-den fractions
%           R = P.*Q   or   R = TIMES(P,Q)
%
% For scalar-den fractions P,Q the command returns scalar-den fraction R
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
% The result is standardly of class 'sdf'. However, the class of the
% result may be specified by an optional input argument CLASS, with
% values 'rdf', 'ldf', 'mdf' or 'sdf'.
%
% See also SDF/MTIMES.

%       Author:  J. Jezek  26-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 26-Apr-2000 $
%                     $ Date 29-May-2000 $
%                     $ Date 06-Nov-2000 $
%                     $ Date 24-Jan-2002 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

tol = PGLOBAL.ZEROING; type = 'sdf';
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
eval('PP = sdf(PP,tol);','error(''Invalid 1st argument.'');');
eval('QQ = sdf(QQ,tol);','error(''Invalid 2nd argument.'');');

[tv,Rv,PP,QQ] = testvf(PP,QQ);
if ~tv, warning('Inconsistent variables.');
end;
[th,Rh,PP,QQ] = testhf(PP,QQ,Rv);
if ~th, warning('Inconsistent sampling periods.');
end;
[td,PP,QQ] = testdf(PP,QQ);
if ~td, error('Matrices not of the same dimensions.');
end;

Rn = 0; Rd = 0;
eval(['Rn = times(PP.frac.num,QQ.frac.num,tol);',...
      'Rd = times(PP.frac.den,QQ.frac.den,tol);'], ...
   'error(peel(lasterr));');
R = sdf(Rn,Rd);
PGLOBAL.COPRIME = cop;

if strcmp(PP.frac.p,'prop') & strcmp(QQ.frac.p,'prop') & ...
      PP.frac.tp==QQ.frac.tp,
   props(R,'prop',PP.frac.tp);
end;

if strcmp(type,'mdf'), R = mdf(R,tol);
elseif strcmp(type,'rdf'), R = rdf(R,tol);
elseif strcmp(type,'ldf'), R = ldf(R,tol);
elseif ~strcmp(type,'sdf'),
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

%end .. @sdf/times

