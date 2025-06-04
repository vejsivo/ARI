function [Y,X,D,E] = stab(B,A,arg3,arg4);
%STAB  Stabilizing controllers for a linear system
%
% The commmands
%    [NC,DC] = STAB(N,D)
%    [NC,DC] = STAB(N,D,'l') 
% compute the right polynomial matrix fraction description 
%    K = -NC/DC 
% of a stabilizing controller for the linear system given by the left 
% coprime polynomial matrix fraction P = D\N. 
%
% The commmand
%    [NC,DC] = STAB(N,D,'r') 
% computes the polynomial matrix fraction description K = -DC\NC of 
% a stabilizing controller of the system given by the right coprime
% fraction P = N/D.
%
% The commmand
%  [NC,DC,E,F] = STAB(N,D), 
% possibly combined with the 'l' or 'r' option, additionally returns 
% the parametrization of all stabilizing controllers. Any stabilizing 
% controller is a minimal realization of the polynomial matrix 
% fraction 
%     K = -Y/X = -(NC*P+E*T)/(DC*P-F*T) 
% or
%	   K = -X\Y = -(P*DC-T*F)\(P*NC+T*E)
% respectively, where T is an arbitrary polynomial matrix and P is
% a stable polynomial matrix such that Y/X (or X\Y, respectively) is 
% proper.
%
% A tolerance TOL may be specified as an additional input argument.
%  
% See also PPLACE, DEBE.

%	Author(s): M. Hromcik, M. Sebek 14-10-98
%	Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 1.0 $  $Date: 14-Oct-1998 10:28:34   $

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');
tol = PGLOBAL.ZEROING;
opt = 'l';

switch nargin
   
  case {0,1},
    error('Not enough input arguments.');
  case 2,
    
  case 3,
    if isnumeric(arg3) 
     tol = arg3;
    elseif ischar(arg3),
     opt = arg3;
    else 
     error('Invalid 3rd argument.');
  end;
  
  case 4,
    if isnumeric(arg3) 
     tol = arg3;
    elseif ischar(arg3),
     opt = arg3;
    else 
     error('Invalid 3rd argument.');
    end;           
    if isnumeric(arg4) 
     tol = arg4;
    elseif ischar(arg4),
     opt = arg4;
    else 
     error('Invalid 4th argument.');
    end;

  otherwise error('Too many input arguments.');
end;	% switch  

if isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;
if ~strcmp(opt,'r') & ~strcmp(opt,'l'),
   error('Invalid command option.');
end;

eval('A = pol(A); B = pol(B);', 'error(peel(lasterr));');

% Check dimensions
eval('[B,A] = testdnd(B,A,opt);', ...
   'error(peel(lasterr));');

% A must be regular:
As = A.s; Bs = B.s;
rankA = rank(A, tol);  
if rankA < As(1),
  error('Denominator matrix is singular.');
end;

% Check variable symbols and sampling periods
var = ''; h = 0;
eval('[var,h,B,A] = testvhnd(B,A);', ...
   'error(peel(lasterr));');

C = pol(zeros(As));
switch opt,
  case 'l',
    % A and B coprime?
    if ~isprime(horzcat(A,B), 'l', tol),
       error('Numerator and denominator are not left coprime.');
    end;
    
    alpha = deg(A, 'row');
    [E,D] = lmf2rmf(B,A);
    
    if As(1)>0,
      delta = deg(D, 'col');
      xi = delta(1) - 1;
	
      gama = alpha + xi;
      gama(gama<0) = 0;

      for i=1:As(1),
        Ci = prand(gama(i), 1, 'sta', var);
        C(i,i) = Ci;
      end;
    end;
   
    eval('[Y,X] = pplace(B,A,C,''l'',tol);', ...
      'error(peel(lasterr));');
     
    %Yc = Y.c;
    %Yc = Yc(:);
    %if all(~Yc(:)), 
    %  Y = pol(zeros(Y.s));
    %  X = pol(eye(X.s));
    %  warning(sprintf(['The system itself is stable. No stabilizing controller is needed.\n', ...
    %  		' 	 To modify the system''s dynamics see PPLACE, DEBE.']));
    %end;  		
    
  case 'r',
    % A and B coprime?
    if ~isprime(vertcat(A,B), 'r', tol),
       error('Numerator and denominator are not right coprime.');
    end;
  
    alpha = deg(A, 'col');
    [E,D] = rmf2lmf(B,A);
    
    if As(1)>0,
      delta = deg(D, 'row');
      xi = delta(1) - 1;
	
      gama = alpha + xi;
      gama(gama<0) = 0;
  
      for i=1:As(1)
        Ci = prand(gama(i), 1, 'sta', var);
        C(i,i) = Ci;
      end;
    end;
   
    eval('[Y,X] = pplace(B,A,C,''r'',tol);', ...
      'error(peel(lasterr));');
    
    %Yc = Y.c;
    %Yc = Yc(:);
    %if all(~Yc(:)), 
    %  Y = pol(zeros(Y.s));
    %  X = pol(eye(X.s));
    %  warning(sprintf(['The system itself is stable. No stabilizing controller is needed.\n', ...
    %  		'	 To modify the system''s dynamics see PPLACE, DEBE.']));
    %end;  	
  
  otherwise error('Invalid command option.');

end;    

%end .. stab
