function P = mat2pol(M)
%MAT2POL  Convert from MATLAB format to Polynomial Toolbox format
%
% The command
%   P = MAT2POL(M)
% converts the polynomial matrix M in MATLAB or Control System Toolbox
% format to a polynomial matrix P in Polynomial Toolbox format.
%
% In MATLAB format a polynomial is represented as a row vector whose
% entries are the coefficients of the polynomial according to descending
% powers. In Control System Toolbox format a polynomial matrix is a cell array
% of the same dimensions as the polynomial matrix, with each cell a
% row vector representing the corresponding entry of the polynomial
% matrix in MATLAB format. If the polynomial matrix consists of a single
% element then in Control System Toolbox format it is represented as a row vector
% with entries in MATLAB format. To comply with an older standard in the
% Control Toolbox, polynomial matrices consisting of a single column may
% be represented as a matrix whose rows are the polynomial entries in
% MATLAB format. The macro POL2MAT converts to this older format if the
% legacy option 'leg' is used. The macro MAT2POL automatically converts
% this older format to a polynomial matrix object.
%
% See also MAT2POL, POL2DSP, DSP2POL.

%     Authors: H. Kwakernaak, R.C.W. Strijbos, November 13, 1998
%     Copyright by PolyX Ltd, 1998
%     $ Revision $  $ Date 22-Jul-2001  J.Jezek  arg check  $
%                   $ Date 25-Jul-2002  J.Jezek  arg check  $

if nargin<1,
   error('Not enough input arguments.');
end;
if isnumeric(M) & ndims(M)==2,
   if isempty(M),
      P = pol(M);
   else
      P = pol(fliplr(M),size(M,2)-1);
   end;
elseif iscell(M)
   P = pol(zeros(size(M)));
   for i = 1:size(M,1)
      for j = 1:size(M,2)
         p = M{i,j};
         P(i,j) = pol(fliplr(p),size(p,2)-1);
      end  
   end
else
   error('Invalid argument; not a polynomial or polynomial matrix in MATLAB format.')
end

%end .. mat2pol

