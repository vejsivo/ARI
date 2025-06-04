function Q = inv(M,arg2,arg3)
%INV    Inverse of scalar-den fraction
%
% The command  Q = INV(M)  returns the inverse of square
% scalar-den fraction P. The result is also a scalar-den fraction.
%
% An optional input argument TOL may specify a zeroing tolerance
% to be used instead the standard one.
%
% An optional input argument MET may specify the method used.
% The methods are  'def','int',  see POL/ADJ.

%       Author:  J. Jezek  26-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 26-Apr-2000 $
%                     $ Date 06-Nov-2000 $
%                     $ Date 06-Jul-2001 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

tol = PGLOBAL.ZEROING;
met = 'int';

if nargin>=2,
   if isa(arg2,'char'), met = arg2;
   elseif isa(arg2,'double'), tol = arg2;
   else error('Invalid 2nd argument.');
   end;
end;
if nargin==3,
   if isa(arg3,'char'), met = arg3;
   elseif isa(arg3,'double'), tol = arg3;
   else error('Invalid 3rd argument.');
   end;
end;

if ~isempty(tol),
   if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;   

Adj = 0; Det = 0; Q = 0;
Mh = M.frac.h;
eval('[Adj,Det] = adj(M.frac.num,tol,met);', ...
   'error(peel(lasterr));');
if eq(Det,0,tol),
   error('Matrix is singular.');
end;
Num = mtimes(M.frac.den,Adj,tol);
Q = sdf(Num,Det);
Q.frac.h = Mh;

if strcmp(PGLOBAL.COPRIME,'cop'),
   Q = coprime(Q,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   Q = reduce(Q);
else Q = smreduce(Q);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'),
   Q = defract(Q);
end;

%end .. @sdf/inv
