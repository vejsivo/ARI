function P = dsp2pol(M)
%DSP2POL Convert from DSP format to polynomial toolbox format
%
% Syntax
%   P = DSP2POL(M)
% In DSP format a polynomial is represented as a row vector whose
% entries are the coefficients of the polynomial according to ascending
% powers. A polynomial matrix is represented as a cell array of the same 
% dimensions as the polynomial matrix, with each cell a row vector
% representing the corresponding entry of the polynomial matrix in 
% DSP format.
%
% See also POL2DSP, POL2MAT, MAT2POL.

%     Author: H. Kwakernaak, September, 1998
%     Copyright by PolyX Ltd, 1998
%     $ Revision $  $ Date 22-Jul-2001  J.Jezek  arg check  $

if nargin<1,
   error('Not enough input arguments.');
end;
if isnumeric(M) & ndims(M)==2 & isempty(M)
   P = pol(M);
elseif isnumeric(M) & ndims(M)==2 & size(M,1)==1 
   P = pol(M,size(M,2)-1);
elseif iscell(M)
   P = pol(zeros(size(M)));
   for i = 1:size(M,1)
      for j = 1:size(M,2)
         p = M{i,j};
         P(i,j) = pol(p,size(p,2)-1);
      end  
   end
else
   error('Invalid argument; not a polynomial or polynomial matrix in DSP format.')
end

%end .. dsp2pol
