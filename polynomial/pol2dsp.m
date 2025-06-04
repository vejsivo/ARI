function M = pol2dsp(P)
%POL2DSP  Convert a polynomial or polynomial matrix to DSP format
%
% Syntax:
%   M = POL2DSP(P)
% In DSP format a polynomial is represented as a row vector whose
% entries are the coefficients of the polynomial according to ascending
% powers. A polynomial matrix is represented as a cell array of the same 
% dimensions as the polynomial matrix, with each cell a row vector
% representing the corresponding entry of the polynomial matrix in 
% DSP format.
%
% See also DSP2POL, POL2MAT, MAT2POL.

%    Author: H. Kwakernaak, September, 1998
%    Copyright by PolyX Ltd, 1998
%    $ Revision $  $ Date 22-Jul-2001  J.Jezek   arg check  $

if nargin<1,
   error('Not enough input arguments.');
end;
eval('P = pol(P);', 'error(peel(lasterr));');

if isempty(P)
   M = zeros(size(P)); return
elseif deg(P) < 0
   if size(p,1) == 1 & size(P,2) == 1
      M = 0;
   else
      M = cell(size(P,1),size(P,2));
      for i = 1:size(P,1)
         for j = 1:size(P,2)
            M{i,j} = 0;
         end
      end
   end
elseif size(P,1) == 1 & size(P,2) == 1
   M  = P{0:deg(P)};
else
   M = cell(size(P,1),size(P,2));
   Pc = P.c;
   p = zeros(1,deg(P)+1);
   for i = 1:size(P,1)
      for j = 1:size(P,2)
         for k = 1:deg(P)+1
            p(1,k) = Pc(i,j,k); 
         end
         M{i,j} = p;
      end
   end
end

%end .. pol2dsp
