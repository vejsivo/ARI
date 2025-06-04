function [a,b,c,d] = abcd(varargin)
%ABCD   State space realization
%
% For expression F convertible to fraction, the command
%    [a,b,c,d] = abcd(F)
% returns the (generalized) observer form realization a,b,c,d,
% that is   F(x) = c*(xI-a)^-1*b + d(x)
% where x is 's','p' (continuous time) or 'z','q' (discrete time)
% or        F(x) = c*x*(I-a*x)^-1*b + d(x^-1)
% where x is 'z^-1','d' (discrete time).
%
% The command
%    [a,b,c,d] = abcd(A,B,C,D,E)
% converts a descriptor system to a (genealized) state space
% system. The systems are equivalent in the sense
%    C*(xE-A)^-1*B + D = c*(xI-a)^-1*b + d(x)
% where x is 's','p' (continuous time) or 'z','q' (discrete time)
% depending on global properties setting.
% Matrix E or matrices D,E may be omitted, the default for E
% being eye matrix, for D zero matrix of corresponding dimensions.
%
% The resulting a,b,c are numerical matrices, d is a numerical
% matrix or a polynomial matrix in 's', 'p', 'z' or 'q'.

% A tolerance TOL may be specified as an additional input argument,
% e.g. 1e-4 . In case of A,B,C... input arguments, the TOL argument,
% to be distinguishable from othr ones, must have a form of a
% character string, e.g. '1e-4' .
% The default value of TOL is the global zeroing tolerance.
%
% See also RDF/ABCD, LDF/ABCD, MDF/ABCD, SDF/ABCD.

%      Author: J. Jezek, 18-Jan-2002
%      Copyright(c) 2002 by Polyx, Ltd.

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
   if isa(F,'ss'),  % LTI argument
      if isempty(F.e),
         a = F.a; b = F.b; c = F.c; d = F.d;
      else
         a = F.a; b = F.b; c = F.c; d = F.d; e = F.e;
         [a,b,c,d] = dss2ss(a,b,c,d,e,tol);
      end;
   else             % fraction argument
      eval('F = sdf(F);', ...
         'error(''Argument is not convertible to fraction.'');');
      eval('[a,b,c,d] = abcd(F,tol);', ...
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
      if size(A,2)~=n | size(B,1)~=n | size(C,2)~=n | ...
            size(D,1)~=p | size(D,2)~=m,
         error('Matrices of inconsistent dimensions.');
      end;
      a = A; b = B; c = C; d = D;
      
   elseif ni==5,
      D = varargin{4}; E = varargin{5};
      if ~isa(D,'double') | ~isa(E,'double'),
         error('Invalid 4th or 5th argument.');
      end;
      eval('[a,b,c,d] = dss2ss(A,B,C,D,E,tol);', ...
         'error(peel(lasterr));');
      
   else
      error('Too many input arguments.');
   end;
end;

%end .. abcd

      
   
   