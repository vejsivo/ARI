function An = power(A,n,arg3,arg4);
%POWER (.^)   Element-wise power of scalar-den fraction
%           A.^n
%
% AN = A.^N or AN = POWER(A,N) denotes element-by-element powers. 
% The scalar-den fraction A and the integer matrix N must have the
% same sizes unless one is a scalar. The scalar operates with
% every element of the other matrix.
%
% If any element of matrix N is negative then the corresponding
% element of A must be nonzero.
%
% AN = POWER(A,N,TOL) works with zeroing specified by the input
% relative tolerance TOL.
%
% The result is standardly of class 'sdf'. However, for N<0, it is
% less time consuming to reach the result in class 'mdf'. For that,
% an optional input argument CLASS may be used. It must be
% a character string, possible values 'mdf', 'sdf', 'rdf' or 'ldf'.

% See also FRAC/MPOWER, SDF/TIMES.

%       Author:  J. Jezek  27-Jan-2000
%       Copyright(c) by Polyx, Ltd.
%       $ Revision $  $ Date 26-Apr-2000 $
%                     $ Date 06-Nov-2000 $
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

if ~isa(n,'double') | ndims(n)>2 | ~isreal(n) |  ...
      (~isempty(n) & any(any(floor(n)~=n)) ),
   error('Invalid power; must be integer.');
end;

An = 0;
if any(size(n)~=1) | n<0,
   eval('An = power(mdf(A),n,tol);','error(peel(lasterr));');
   if strcmp(type,'mdf'),
   elseif strcmp(type,'sdf'), An = sdf(An,tol);
   elseif strcmp(type,'rdf'), An = rdf(An,tol);
   elseif strcmp(type,'ldf'), An = ldf(An,tol);
   else error('Invalid command option.');
   end;
else
   eval('An = sdf(power(A.frac.num,n,tol),power(A.frac.den,n,tol));', ...
      'error(peel(lasterr));');
   if strcmp(type,'sdf'),
   elseif strcmp(type,'mdf'), An = mdf(An,tol);
   elseif strcmp(type,'rdf'), An = rdf(An,tol);
   elseif strcmp(type,'ldf'), An = ldf(An,tol);
   else error('Invalid command option.');
   end;
end;

allnpos = all(all(n>=0));
if strcmp(A.frac.p,'prop') & allnpos,
   props(An,'prop',A.frac.tp);
end;
if strcmp(A.frac.c,'cop'),
   props(An,'cop',A.frac.tc);
end;
if strcmp(A.frac.r,'red') & allnpos,
   props(An,'red');
end;

if strcmp(PGLOBAL.COPRIME,'cop'),
   An = coprime(An,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   An = reduce(An);
else
   An = smreduce(An);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'),
   An = defract(An);
end;

%end .. @sdf/power
