function [N,D] = ss2rmf(varargin)
%SS2RMF  Convert a state space system to a right polynomial matrix fraction
%
% The commands
%    [N,D] = SS2RMF(a,b,c)
%    [N,D] = SS2RMF(a,b,c,d)
%    [N,D] = SS2RMF(a,b,c,d,tol)
% with d a constant matrix or a polynomial matrix, calculate a
% right coprime polynomial fraction representation
%     N D^-1
% such that
%     N D^-1 = c (xI-a)^-1 b + d(x),
% where x is 's','p' (continuous-time) or 'z','q' (discrete-time)
% or
%     N D^-1 = c x (I - x a)^-1 b + d(x^-1),
% where x is 'z^-1','d' (discrete time),
% and the denominator matrix D is column reduced.
% The parameter tol is a tolerance that is used in determining 
% the controllability indexes (that is, the column degrees of D).

% Author: Rens C.W. Strijbos, December 08, 1998.
% Copyright 1998 by Polyx, Ltd.
% $ Revision 3.0 $  $ Date 12-Oct-2000  J.Jezek $ 
%                   $ Date 28-Jul-2001  J.Jezek $
%                   $ Date 24-Feb-2003  J.Jezek $

switch nargin
case {3,4,5}
   a = varargin{1}; b = varargin{2}; c = varargin{3};
   if ~isa(a,'double') | ~isa(b,'double') | ~isa(c,'double') | ...
         ndims(a)>2 | ndims(b)>2 | ndims(c)>2,
      error('Invalid 1st, 2nd or 3rd argument; must be const.');
   end;
   [ra,ca] = size(a); [rb,cb] = size(b); [rc,cc] = size(c);
   switch nargin
   case 3
      d = pol(zeros(rc,cb));
      tol = eps*10^2;
   case {4,5}
      d = varargin{4};
      if isa(d,'double') & ndims(d)==2,
         d = pol(d);
      elseif ~isa(d,'pol'),
         error('Invalid 4th argument; must be const or pol.');
      end;
      switch nargin
      case 4
         tol = eps*10^2;
      case 5
         tol = varargin{5};
         if ~(isa(tol,'double') & ~any(size(tol)-[1 1]) & ...
               isreal(tol) & tol>=0 & tol<=1),
            error('Invalid tolerance.')
         end
      end
   end
case {1,2}
   error('Not enough input arguments.')
otherwise
   error('Too many input arguments.');
end

var = d.var;
if ~isempty(var),
   if ~(strcmp(var,'s') | strcmp(var,'p') | ...
        strcmp(var,'z') | strcmp(var,'q')),
      error('Invalid variable symbol in 4th argument.');
   end;
end;  

[rd,cd] = size(d);
if (ra~=ca | ra~=rb | ra~=cc | cb~=cd | rc~=rd),
   error('Matrices of incompatible dimensions.');
end;

[f,g,h,n] = bhf(a',c',b',tol);
[F,G,H,n] = bhf(f',h',g',tol);

if length(n) == 0
   D = pol(eye(cb),0);
   N = pol(zeros(rc,cb),0);
else
   [N,D] = bhf2rmf(F,G,H,n,tol);
end;

if ~isempty(var),
   N.var = var;
   D.var = var;
end
N = N + d * D;
%if strcmp(var, 'z^-1') | strcmp(var, 'd'),
%   [N,D] = reverse(N,D,'r',tol);
%end;
%  commented out by J.Jezek 24-Feb-2003
%  as it can never happen, see above var=d.var
%  see also corrections in SS2LMF, BHF2RMF

%end .. ss2rmf

