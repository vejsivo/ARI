function M = pol2mat(P,arg)
%POL2MAT  Convert from Polynomial Toolbox format to MATLAB format
%
% Syntax:
%    M = POL2MAT(P[,'LEG'])
% In MATLAB format a polynomial is represented as a row vector whose
% entries are the coefficients of the polynomial according to descending
% powers. In Control System Toolbox format a polynomial matrix is a cell array
% of the same dimensions as the polynomial matrix, with each cell a
% row vector representing the corresponding entry of the polynomial
% matrix in MATLAB format. If the polynomial matrix consists of a single
% element then in Control System Toolbox format it is represented as a row vector
% with entries in MATLAB format. To comply with an older standard in the
% Control System Toolbox, polynomial matrices consisting of a single column may
% be represented as a matrix whose rows are the polynomial entries in
% MATLAB format. The macro POL2MAT converts to this older format if the
% legacy option 'leg' is used. The macro MAT2POL automatically converts
% this older format to a polynomial matrix object.
%
% See also MAT2POL, POL2DSP, DSP2POL.

%    Authors: H. Kwakernaak, R.C.W. Strijbos, December 3, 1998
%    Copyright by PolyX Ltd, 1998
%    $ Revivion $  $ Date 22-Jul-2001  J.Jezek  arg check  $

if nargin < 1
   error('Not enough input arguments.');
end;
legacy = 0;
if nargin > 1
   if ischar(arg) & strcmp(arg,'leg')
      legacy = 1;
   else
      error('Invalid command option.')
   end
end

eval('P = pol(P);', 'error(peel(lasterr));');
if isempty(P)
   M = zeros(size(P)); return
elseif deg(P) < 0
   if size(P,2) == 1 & legacy
      M = zeros(size(P));
   else
      M = cell(size(P,1),size(P,2));
      for i = 1:size(P,1)
         for j = 1:size(P,2)
            M{i,j} = 0;
         end
      end
   end
elseif length(P) == 1 | (size(P,2) == 1 & legacy)
   M  = P{deg(P):-1:0};
else
   M = cell(size(P,1),size(P,2));
   Pc = P.c;
   p = zeros(1,deg(P)+1);
   for i = 1:size(P,1)
      for j = 1:size(P,2)
         for k = 1:deg(P)+1
            p(1,k) = Pc(i,j,k); 
         end
         M{i,j} = fliplr(p);
      end
   end
end

%end .. pol2mat
