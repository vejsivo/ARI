function [a,b,c,d,e] = abcde(varargin)
%ABCDE   Descriptor state space realization
%
% For expression F convertible to fraction, the command
%    [a,b,c,d,e] = ABCDE(F)
% returns the descriptor state space realization a,b,c,d,e, that is
% that is   F(x) = c*(xe-a)^-1*b + d
% where x is 's','p' (continuous time) or 'z','q' (discrete time)
% or        F(x) = c*x*(e-a*x)^-1*b + d
% where x is 'z^-1','d' (discrete time).
%
% The command
%    [a,b,c,d,e] = abcde(A,B,C,D)
% converts a (generalized) state space system to a descriptor
% system. The systems are equivalent in the sense
%     C*(xI-A)^-1*B + D(x) = c*(xe-a)^-1*b + d
% where x is 's','p' (continuous time) or 'z','q' (discrete time).
% Matrix D may be onitted, the default being the zero matrix of
% corresponding dimensions.
%
% The resulting a,b,c,d,e are numerical matrices.
%
% A tolerance TOL may be specified as an additional input argument,
% e.g. 1e-4. In case of A,B,C... input arguments, the TOL argument,
% to be distinguishable from other ones, must have a for of a
% character string, e.g. '1e-4' .
% The default value of TOL is the global zeroing tolerance.
%
% See also RDF/ABCDE, LDF/ABCDE, MDF/ABCDE, SDF/ABCDE.

%     Author: J. Jezek, 10-Mar-2000
%     Copyright (c) 2000 by Polyx, Ltd.
%     $ Revision $  $ Date 15-Jun-2000 $
%                   $ Date 20-Jan-2002 $

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
   
elseif ni<=2,
   if ni==1,
      tol = PGLOBAL.ZEROING;
   else
      tol = varargin{2};
      if ~isempty(tol) & ~isa(tol,'double'),
         error('Invalid tolerance.');
      end;
   end;
   
   F = varargin{1};
   if isa(F,'ss'),     % LTI argument
      a = F.a; b = F.b; c = F.c; d = F.d;
      if isempty(F.e),
         e = eye(size(F.a));
      else
         e = F.e;
      end;
   else                % fraction argument
      eval('F = sdf(F);', ...
         'error(''Argument is not convertible to fraction.'');');
      eval('[a,b,c,d,e] = abcde(F,tol);', ...
         'error(peel(lasterr));');
   end;
   
else                % state space arguments
   A = varargin{1}; B = varargin{2}; C = varargin{3};
   if ~isa(A,'double') | ~isa(B,'double') | ~isa(C,'double'),
      error('Invalid 1st, 2nd or 3rd argument.');
   end;
   n = size(A,1); m = size(B,2); p = size(C,1);
   
   tol = PGLOBAL.ZEROING;
   if ni>3,
      lastarg = varargin{ni};
      if isa(lastarg,'char'),
         if ~isempty(lastarg),
            tol = str2num(lastarg);
            if isempty(tol),
               error('Invalid last argument.');
            end;
            if ~(length(tol)==1 & isreal(tol) & tol>=0 & tol<=1),
               error('Invalid tolerance.');
            end;
         end;
         ni = ni-1;
      elseif ni==6 & isa(lastarg,'double'),
         if ~isempty(lastarg),
            tol = lastarg;
            if ~(length(tol)==1 & isreal(tol) & tol>=0 & tol<=1),
               error('Invalid tolerance.');
            end;
         end;
         ni = ni-1;
      end;
   end;
   
   if ni<=4,
      if ni==3,
         D = zeros(p,m);
      else
         D = varargin{4};
         if ~isa(D,'double') & ~isa(D,'pol'),
            error('Invalid 4th argument.');
         end;
         if isa(D,'pol'),
            if strcmp(D.v,'z^-1') | strcmp(D.v,'d'),
               error('Invalid variable symbol in 4th argument.');
            end;
         end;
      end;
      eval('[a,b,c,d,e] = ss2dss(A,B,C,D,tol);', ...
         'error(peel(lasterr));');
      
   elseif ni==5,
      D = varargin{4}; E = varargin{5};
      if ~isa(D,'double') | ~isa(E,'double'),
         error('Invalid 4th or 5th argument.');
      end;
      if size(A,2)~=n | size(B,1)~=n | size(C,2)~=n | ...
            size(D,1)~=p | size(D,2)~=m | size(E,1)~=n | ...
            size(E,1)~=n,
         error('Matrices of inconsistent dimensions.');
      end;
      a = A; b = B; c = C; d = D; e = E;
      
   else
      error('Too many input arguments.');
   end;
end;

%end .. abcde




