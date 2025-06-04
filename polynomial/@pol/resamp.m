function B = resamp(A,k,h)
%RESAMP     Resample discrete-time polynomial
%
% Let A be a polynomial in variable 'z^-1' or 'd'
%    A = A(0) + ... + A(n)*z^-n + ...
% whose n-th coefficient is A(n). The command
%    B = RESAMP(A,K,H)
% where scalar integers K,H are resampling period, K>0, and
% resampling phase, H>=0, (defaults being K=2, H=0), returns
% polynomial B in the same variable, with n-th coefficient
%  B(n) = A(K*n+H). Note that the meaning of the variable
% has changed:  new var = old var to the K-th power.
% The sampling period  B.h  of B is K-times greater than  A.h . 
%
% For example:
%       A = 1 + 0.8*z^-1 + 0.64*z^-2 + 0.512*z^-3 +
%           0.4096*z^-4 + 0.32768*z^-5
%       RESAMP(A,3,0) = 1 + 0.512*z^-1
%       RESAMP(A,3,1) = 0.8 + 0.4096*z^-1
%       RESAMP(A,3,2) = 0.64 + 0.32768*z^-1
%
% For a polynomial in variable 'z' or 'q', the powers of the variable
% grow in the opposite direction to the numbering of H. The formula
% is  B(n) = A(K*n-H). This rule is for compatibility with TSP/RESAMP.
%
% For example:
%     A = 0.32768*z^5 + 0.4096*z^4 + 0.512*z^3 +
%         0.64*z^2 + 0.8*z + 1
%     RESAMP(A,3,0) = 0.512*z + 1
%     RESAMP(A,3,1) = 0.32768*z^2 + 0.64*z
%     RESAMP(A,3,2) = 0.4096*z^2 + 0.8*z
%
% For a polynomial in variable 's' or 'p', the resampling is
% not defined.
%
% H can also be, instead of scalar, a vector. In such
% a case, the result is a cell vector.
%
% See also TSP/RESAMP, RDF/RESAMP, LDF/RESAMP, POL/DILATE.

%        Author:  J. Jezek 08-Nov-1999
%        Copyright(c) 1999 by Polyx, Ltd.
%        $ Revision $  $ Date 19-Jun-2000 $
%                      $ Date 04-Oct-2000 $
%                      $ Date 01-Feb-2001 $
%                      $ Date 28-Feb-2003 $
                      
ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ~isa(A,'pol'),
   error('Some argument but not 1st is invalidly polynomial.');
end;

if ni<2 | isempty(k),
   k = 2;
else
   if ~isa(k,'double') | length(k)~=1 | ~isreal(k) | ...
         floor(k)~=k | k<=0,
      error('Invalid resampling period.');
   end;
end;

if ni<3 | isempty(h),
   h = 0;
else
   if isa(h,'double') & ndims(h)==2 & any(size(h)==1) & ...
         isreal(h) & all(floor(h)==h) & all(h>=0),
   else error('Invalid resampling phase.');
   end;
end;

Av = A.v;
if strcmp(Av,'s') | strcmp(Av,'p'),
   error('Invalid variable symbol; must be discrete-time.');
elseif strcmp(Av,'z') | strcmp(Av,'q'),
   h = -h;
end;

lh = length(h);
if lh>1, B = cell(size(h));
end;

Ad = A.d; As = A.s; Ah = A.h;
for i = 1:lh,
   hi = h(i);
   BB.d = floor((Ad-hi)/k);
   if BB.d<0, BB.d = -Inf;
   end;
   BB.s = As;
   if hi>=0,
      BB.c = A.c(:,:,1+hi:k:end);
   else
      BB.c = cat(3,zeros(As(1),As(2)),A.c(:,:,1+k+hi:k:end));
   end;
   BB.v = Av;
   if isempty(Ah), BB.h = Ah;
   else BB.h = k*Ah;
   end;
   BB.u = [];
   BB.version = 3.0;

   BB = class(BB,'pol');
   BB = pclear(BB);
   
   if lh==1, B = BB;
   else B{i} = BB;
   end;
   clear BB;
end;  %for i

%end .. @pol/resamp

