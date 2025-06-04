function P = root2pol(Z, K, arg3, arg4);
%ROOT2POL  Construct a polynomial matrix from its zeros and gains.
%
% Given a cell array Z of vectors and a matrix K the command
%    P = ROOT2POL(Z,K)
% returns the polynomial matrix P whose entries P(i,j) have zeros Z{i,j} 
% and gains K(i,j). If K is missing then it is supposed to be a ones-array. 
% The commands
%    P = ROOT2POL(Z,K,VAR)
%    P = ROOT2POL(Z,VAR) 
% allow the user to specify the variable of P.
%    
% A tolerance TOL may be specified as an additional input argument.
%
% See also POL, ROOTS.

%	Author(s): M. Hromcik, M. Sebek 15-9-98
%	Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 1.0 $  $Date: 15-Sep-1998 10:28:34   $
%       $Revision: 2.0 $  $Date: 10-Oct-2000 10:30:00 - Version 3.0, M. Hromcik $
%                         $Date: 10-Aug-2001  J.Jezek   arg checking     $

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;'); 
tol = PGLOBAL.ZEROING;

Var = [];
if nargin==0,
   error('Not enough input arguments.');
end; 
if ~isa(Z, 'cell'),
   error('Invalid 1st arguent; must be a cell array.');
end;
Zs = size(Z);

switch nargin,
  case 1,
    K = ones(Zs);
     
  case 2,  
    if isempty(K), K = ones(Zs); 
    elseif isa(K, 'char'), Var = K; K = ones(Zs);
    elseif ~isnumeric(K),
       error('Invalid 2nd argument.');
    end;
    
  case 3,
    if isempty(K), K = ones(Zs);
    elseif ~isnumeric(K),
       error('Invalid 2nd argument.');
    end;
    if isa(arg3,'char'), Var = arg3;
    elseif isnumeric(arg3), tol = arg3;
    else error('Invalid 3rd argument.');
    end;
    
  case 4,
    if isempty(K), K = ones(Zs);
    elseif ~isnumeric(K),
       error('Invalid 2nd argument.');
    end;    
    if isa(arg3,'char'), Var = arg3;
    elseif isnumeric(arg3), tol = arg3;
    else error('Invalid 3rd argument.');
    end;
    if isa(arg4,'char'), Var = arg4;
    elseif isnumeric(arg4), tol = arg4;
    else error('Invalid 4th argument.');
    end;
end;
 
Ks = size(K);
if any(Zs-Ks),
   error('1st and 2nd argument must have the same dimensions.');
end;

P = pol(zeros(Zs));
if Zs(1)>0 & Zs(2)>0,
  for i = 1:Zs(1),
    for j = 1:Zs(2),
       Zij = Z{i,j};
       if ~isnumeric(Zij) | ndims(Zij)>2 | all(size(Zij)>1), 
         error('Invalid entry in 1st argunent; must be a numeric vector.'); 
       end;
       P(i,j) = K(i,j) * pol( fliplr(poly(Zij)), length(Zij) );
    end; 
  end;   
end;

if ~isempty(Var),
   eval('P.v = Var;', 'error(peel(lasterr));');
end;

% Zeroing:
% Pc = P.c;
% me = max(abs(Pc(:)));
% P = pzer(P, me*tol);        ???  Now the tolerance is not used at all

%end .. root2pol
