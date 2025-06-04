function F = frac(N,D)
%FRAC     Create fraction
%
% The command  F = FRAC(N,D) creates fraction F from
% given polynomial matrix N and polynomial matrix D.
% Polynomials N,D must be finite (Inf, NaN not allowed).
%
% The command  F = FRAC(N) returns N if N is aleady a fraction.
% Otherwise, it means  F = FRAC(N,1).
%
% The variable symbols of polynomials N,D should be the same. When they
% are not, a warning is issued and the symbols are changed to
% PGLOBAL.VARIABLE. However, if one symbol is 'z' and the other 'z^-1'
% then the symbols play a role, no warning being issued, the resulting
% symbol is taken from N.
%
% This macro is called only within the constructors of
% objects "left-denominator fraction", "right-denominator
% fraction", "matrix-denominator fraction" or "scalar-
% denominator fraction". The class "fraction" is auxiliary
% and its purpose is only to serve as a common ancestor
% of the above mentioned classes. No objects of class
% "fraction" should be created.

% The structure of object "fraction" is as follows:
%
% M.num      numerator, polynomial matrix
% M.den      denominator, polynomial matrix
% M.s        size
% M.v        variable symbol
% M.h        sampling period
% M.u        user data
% M.c        flag of 'coprime': 'cop', 'ncop' or 'cop?'
% M.r        flag of 'reduced': 'red', 'nred' or 'red?'
% M.p        flag of 'proper' : 'prop', 'nprop' or 'prop?'
% M.tc       tolerance used for 'cop'
% M.tp       tolerance used for 'prop'
% M.version  version number, always 3.0

%       Author:  J. Jezek,  18-Feb-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 28-Nov-2000 $
%                     $ Date 02-Aug-2001 $
%                     $ Date 22-Jul-2002 $
%                     $ Date 06-Oct-2002 $
%                     $ Date 14-Oct-2002 $

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
Nc = class(N);
if ni==1,
   if strcmp(Nc,'frac'),  % quick return
      F = N; return;
   end;
   if isa(N,'tsp'),
      N = sdf(N);
   end;
   if isa(N,'rdf') | isa(N,'ldf') | isa(N,'mdf') | isa(N,'sdf'),
      D = N.den; N = N.num;
      ni = 0;
   end;
end;

if ni>=1,
   eval('N = pol(N);','error(peel(lasterr));');
end;
if ni==2,
   eval('D = pol(D);','error(peel(lasterr));');
elseif ni==1,
   D = pol(1);
end;
   
[tv,v,N,D] = testvp(N,D);
if tv==2,
   var = pol([0 1],1,v);
   factor = var^deg(D);
   N = N*factor; D = pol(D*factor);
elseif tv==0, warning('Inconsistent variables.');
end;

[th,h,N,D] = testhp(N,D,v);
if ~th, warning('Inconsistent sampling periods.');
end;

if ~all(all(isfinite(N))) | ~all(all(isfinite(D))),
   error('Polynomial is not finite.');
end;

superiorto('double','pol','tsp');

F.num = N;
F.den = D;
F.s = N.s;
if any(F.s==0), v = '';
end;
F.v = v;
if isempty(v), h = [];
end;
F.h = h;
F.u = [];
F.c = 'cop?';
F.r = 'red?';
F.p = 'prop?';
F.tc = [];
F.tp = [];
F.version = 3.0;

F = class(F,'frac');

%end .. @frac/frac
