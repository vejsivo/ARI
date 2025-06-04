function [Act,varargout] = ctranspose(A);
%CTRANSPOSE   Conjugate transpose polynomial
%                A'
%  
% The commands
%    AT = A'  
%    AT = CTRANSPOSE(A) 
% return the conjugate transpose of the polynomial matrix A.
%
% If A is a polynomial matrix in variable 's', then 
%    AT(s) = CONJ(A'(-s)) .
% For the variable 'p' the command works similarly.           
%
% If A is a polynomial matrix in variable 'z', then
%    AT(z) = CONJ(A'(z^-1)) .
% The result is a polynomial in 'z^-1'.
% For the variable 'z^-1' the command works similarly,
% the result is a polynomial in 'z'.
%
% If A is a polynomial matrix in variable 'q', then 
%    AT(q) = CONJ(A'(q^-1)*q^n) ,
% where n is the degree of A. The reault is again
% a polynomial in 'q'. The command
%    [AT,N] = CTRANSPOSE(A) 
% returns both AT and the degree of A.  For the variable 'd' 
% the commands work similarly.
%
% See also POL/TRANSPOSE, POL/CONJ.

%     	Author(s): M. Hromcik, M. Sebek 16-2-98
%     	Copyright (c) 1998 by Polyx, Ltd.
%        $Revision: 2.0 $  $Date: 06-Mar-1998 09:15:34   $
%        $Revision: 3.0 $  $Date: 11-Aug-1999 12:00:00  J.Jezek  $
%                          $Date: 12-Apr-2000 12:00:00  J.Jezek  $

% Effect on other properties:
% Act.u: UserData are deleted.
   
no = nargout;
 
Act.d = A.d;
Act.s = fliplr(A.s);
Act.c = conj( permute(A.c,[2,1,3]) );

% variables 's' and 'p'  
if strcmp(A.v,'s') | strcmp(A.v,'p'),
 if no > 1, 
  error('Too many output arguments.');
 end; 
 for i=1:Act.d+1,
  Act.c(:,:,i)=(-1)^(i-1) * Act.c(:,:,i);
 end;
 Act.v = A.v;
 Act.h = A.h;
 Act.u = [];
 Act.version = 3.0;
 Act = class(Act,'pol');
  
% variables 'q' and 'd' 
elseif strcmp(A.v,'q') |  strcmp(A.v,'d'),
 if no > 2, 
  error('Too many output arguments.');
 end;
 Act.c = flipdim(Act.c,3); 
 Act.v = A.v;
 Act.h = A.h;
 Act.u = [];
 Act.version = 3.0;
 Act = class(Act,'pol');
 Act = pclear(Act); 
 varargout{1} = A.d;
 
%variables 'z' and 'z^-1' 
elseif strcmp(A.v,'z') | strcmp(A.v,'z^-1'),
 if no > 1,
  error('Too many output arguments.');
 end;
 if strcmp(A.v,'z'), Act.v = 'z^-1';
 else Act.v = 'z';
 end;
 
 Act.h = A.h;
 Act.u = [];
 Act.version = 3.0;
 Act = class(Act,'pol');
 
%variable ''
else
 if no > 2,
  error('Too many output arguments.')
 end;
 Act.v = '';
 Act.h = [];
 Act.u = [];
 Act.version = 3.0;
 Act = class(Act,'pol');
 varargout{1} = 0;
   
end;	%if

%end .. @pol/ctranspose    

    

 

 
