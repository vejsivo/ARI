function y = norm(A,arg1,arg2,arg3);
%NORM  Norms of polynomial matrices and polynomial matrix fractions
%
% For a polynomial matrix 
%    P = P0 + P1 var + ... + Pd var^d 
% the commmand
%    NORM(P)
% returns the largest singular value of the constant matrix [P0, P1, ... , Pd].
% The commmand
%    NORM(P,'blk')  
% is the same as NORM(P), while
%    NORM(P,'lead') 
% is the largest singular value of the leading coefficient matrix Pd. The
% commmand
%    NORM(P,'max')  
% is the maximum of largest singular values of the coefficient matrices P0,
% P1, .. ,Pd. The commmand
%    NORM(P,'blk',TYPE)  
% returns the TYPE norm (TYPE is 1, 2, inf, 'fro') of the constant matrix 
% [P0, P1, ... , Pd]. The commmand
%    NORM(P,'lead',TYPE) 
% returns the TYPE norm (TYPE is 1, 2, inf, 'fro') of the highest coefficient 
% matrix of P. The commmand
%    NORM(P,'max',TYPE)  
% returns the maximum of TYPE norms of the coefficient matrices P0, P1, .. ,Pd.
%
% For polynomial matrix fractions
%    NORM(N,D)
%    NORM(N,D,Inf)
%    NORM(N,D,'l',Inf) 
% compute the H-infinity norm of the rational matrix D\N, while
%    NORM(N,D,'r')
%    NORM(N,D,'r',Inf) 
% compute the H-infinity norm of the rational matrix N/D. The commmands
%    NORM(N,D,2)
%    NORM(N,D,'l',2)
%    NORM(N,D,'r',2) 
% compute the H2-norm of the rational matrix D\N or N/D, respectively.
%
% See also POL/HINFNORM, POL/H2NORM.

%	Author(s): M. Hromcik, M. Sebek 16-2-98
%	Copyright (c) 1998-2000 by PolyX, Ltd.
%       $Revision: 4.0 $  $Date: 23-Oct-1998 10:28:34   $
%
%		  $ Version 3 - 22-Nov-2000, Martin Hromcik $
%       $             17-Jul-2001, Jan Jezek      $



type1 = 'blk';
type2 = 2;

switch nargin,
 case 0,
  error('Not enough input arguments.');
 
 case 1,
  y = norm( A.c(:,:) );
  return;

 case 2,
  if isnumeric(arg1) & length(arg1) ==1,
    type2 = arg1;
  elseif isa(arg1,'char'), 
    type1 = arg1;
  else 		% fraction norm.
    y = hinfnorm(A,arg1);   
    return;
  end;  
    
 case 3,
  if isnumeric(arg1) & length(arg1) ==1,
    type2 = arg1;
  elseif isa(arg1,'char'), 
    type1 = arg1;
  else 		% fraction norm
    if ischar(arg2),
      y = hinfnorm(A,arg1,arg2);
      return;
    elseif arg2 == 2,
      y = h2norm(A,arg1);  
      return;
    elseif isinf(arg2),
      y = hinfnorm(A,arg1);
      return;
    else error('Invalid command option.');  
    end;  
  end;
     
  if isnumeric(arg2) & length(arg2) ==1,
    type2 = arg2;
  elseif isa(arg2,'char'), 
    type1 = arg2;
  else
    error('Invalid command option.');  
  end; 
  
  case 4,	% fraction norm
    type = Inf;
  
    if ischar(arg2), opt = arg2;
    elseif arg2 == 2, type = 2;
    elseif isinf(arg2), type = Inf;
    else error('Invalid command option.');
    end;
    
    if ischar(arg3), opt = arg3;
    elseif arg3 == 2, type = 2;
    elseif isinf(arg3), type = Inf;
    else error('Invalid command option.');
    end;
  
    if type == 2,
      y = h2norm(A,arg1,opt);
      return;
    elseif isinf(type), 
      y = hinfnorm(A,arg1,opt);
      return;
    else error('Invalid command option.');
    end;      
     
 otherwise error('Too many input arguments.');
end; 	%switch  

Ac = A.c;

switch type1
 case 'max',			 % 'max' norm	
	 y = 0;	
	 for i=1:A.d+1,
	  y = max( y, norm(Ac(:,:,i), type2) );
	 end;

 case 'lead',			 % 'lead' norm
         y = norm( Ac(:,:,max(end,1)), type2 );
	
 case 'blk',			 % 'block' norm
	 y = norm( Ac(:,:),type2 );
	 
 otherwise 
	 error('Invalid command option.')
end;
 
%end .. @pol/norm
