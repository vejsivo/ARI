function R = rdivide(P,Q,arg3,arg4);
%RDIVIDE(./)  Element-wise right divide scalar-den fractions
%           R = P./Q   or   R = RDIVIDE(P,Q)
%
% For scalar-den fractions P,Q, the command returns scalar-den fraction R
% whose (i,j)-th entry is  P(i,j)/Q(i,j) . All entries of Q must be
% nonzero.
%
% Matrices P,Q must have the same dimensions, unless one of them is
% scalar. In such a case, this scalar operates with every entry
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
% The result is standardly of class 'sdf'. However, it is less time
% consuming to reach the result in class 'mdf'. For that, an optional
% input argument CLASS may be used. It must be a character string,
% possible values 'mdf', 'sdf', 'rdf' or 'ldf'.
%
% See also SDF/LDIVIDE, FRAC/MRDIVIDE, FRAC/MLDIVIDE.

%       Author:  J. Jezek  04-Feb-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 26-Apr-2000 $
%                     $ Date 29-May-2000 $
%                     $ Date 03-Aug-2000 $
%                     $ Date 24=Jan-2002 $
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

eval('P = sdf(P);','error(''Invalid 1st argument.'');');
eval('Q = sdf(Q);','error(''Invalid 2nd argument.'');');
if any(any(Q.frac.num==0)),
   error('Denominator is zero.');
end;

[tv,Rv,P,Q] = testvf(P,Q);
if ~tv, warning('Inconsistent variables.');
end;
[th,Rh,P,Q] = testhf(P,Q,Rv);
if ~th, warning('Inconsistent sampling periods.');
end;
[td,P,Q] = testdf(P,Q);
if ~td, error('Matrices not of the same dimensions.');
end;

num = times(P.frac.num,Q.frac.den,tol);
den = times(P.frac.den,Q.frac.num,tol);
R = mdf(num,den);

if strcmp(type,'sdf'),
   R = sdf(R,tol);
elseif strcmp(type,'ldf'),
   R = ldf(R,tol);
elseif strcmp(type,'rdf'),
   R = rdf(R,tol);
elseif ~strcmp(type,'mdf'),
   error('Invalid command option.');
end;

if strcmp(PGLOBAL.COPRIME,'cop'), R = coprime(R,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'), R = reduce(R);
else R = smreduce(R);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), R = defract(R);
end;

%end .. @sdf/rdivide
