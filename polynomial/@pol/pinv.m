function M = pinv(A, arg2, arg3);
%PINV   Pseudoinverse of a full-rank polynomial matrix
%
% The command
%     M = PINV(A) 
% returns the pseudoinverse of a polynomial matrix A, that is, M such that
% the relations
%   A * M * A = A 
%   M * A * M = M
% hold and A*M and M*A are symmetric (i.e., (A*M).' = A*M, (M*A).' = M*A). 
% A must have full rank, i.e., RANK(A) == MIN(SIZE(A)) must hold.
%
% M is a scalar denominator fraction.
%
% The commands
%    PINV(A,'int') or PINV(A,'def') 
% allow the user to specify the method used. For details about the methods 
% see POL/ADJ.
% 
% The command
%    PINV(A,TOL) 
% uses the optional relative tolerance TOL for zeroing and rank testing. 
% The tolerance and method options can be combined together.
%
% See also POL/INV, POL/RANK.

%       Author(s): M. Hromcik, M. Sebek 5-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 16-Oct-1998 10:28:34   $
%       $Revision: 3.0 $  $Date: 12-Oct-2000            Version 3, M. Hromcik $
%                         $Date: 08-Jul-2001    J.Jezek  $

% Operator ' replaced by .' for possible problems with discrete matrices.

global PGLOBAL;

tol = PGLOBAL.ZEROING;
met = 'def';

switch nargin,
  case 1,
    rankA = rank(A);
  case 2,
    if isnumeric(arg2) & length(arg2)==1 & ...
          isreal(arg2) & arg2>=0 & arg2<=1,
      tol = arg2;
	   rankA = rank(A,tol);
    else
    	met = arg2;
    	rankA = rank(A);
    end;
  otherwise
    usrtol = 0;
    if isnumeric(arg2) & length(arg2)==1 & ...
          isreal(arg2) & arg2>=0 & arg2<=1,
	   tol = arg2;
	   usrtol = 1;
    else
    	met = arg2;
    end;
    if isnumeric(arg3) & length(arg3)==1 & ...
          isreal(arg3) & arg3>=0 & arg3<=1,
	   tol = arg3;
	   usrtol = 1;
    else
    	met = arg3;
    end;
    if usrtol, rankA = rank(A,tol);
    else, rankA = rank(A);
    end;
end;	%switch

As = A.s;

if rankA == min(As),
  
	if As(1) == As(2),
      eval('[Nx, dx] = adj(A, met, tol);', ...
         'error(peel(lasterr));')
	
	elseif As(1) < As(2),
      eval('[Nx,dx] = adj(A*A.'', met, tol);', ...
         'error(peel(lasterr));');
	   Nx = A.'*Nx;

	else
      eval('[Nx,dx] = adj(A.''*A, met, tol);', ...
         'error(peel(lasterr));');
	   Nx = Nx*A.';
	end;  

else			% rank deficient	

	% [B,U]  = tri(A,'col');
	% [C,V]  = tri(A,'row');
	% [T1,W] = tri(C,'row');
	% [T2,Z] = tri(B,'col');
	
	error('Matrix does not have full rank.');
	
end;	

eval('M = sdf(Nx,dx);','error(peel(lasterr));');

isproper(M);
if strcmp(PGLOBAL.COPRIME,'cop'),
   M = coprime(M,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   M = reduce(M);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'),
   M = defract(M);
end;

%end .. @pol/pinv
