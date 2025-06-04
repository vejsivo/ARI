function [Y,X,D,E,degT] = debe(B,A,arg3,arg4);
%DEBE  Deadbeat controllers of discrete-time linear system
%
% The commmands
%   [NC,DC] = DEBE(N,D)
%   [NC,DC] = DEBE(N,D,'l') 
% compute the right polynomial matrix fraction description K = -NC/DC 
% of a deadbeat controller of the linear discrete-time strictly proper 
% system given by the left polynomial matrix fraction P = D\N. 
%
% The commmand
%    [NC,DC] = DEBE(N,D,'r') 
% computes the left polynomial matrix fraction description K = -DC\NC 
% of a deadbeat controller of the discrete-time strictly proper system 
% P = N/D. 
%
% The commmand
%    [NC,DC,E,F,DEGT] = DEBE(N,D), 
% possibly combined with the 'l' or 'r' options, additionally returns 
% the parametrization of the deadbeat controllers. Their transfer 
% functions may be expressed as
%    K = -(NC+E*T)/(DC-F*T) 
% or
%	  K = -(DC-T*F)\(NC+T*E),
% respectively. If the variable of A and B is 'z' or 'q' then T is a 
% free polynomial matrix with entry degrees less than or equal to DEGT. 
% For the variables 'd' and 'z^-1', T is a free polynomial matrix
% such that X{0} is nonsingular.
%
% A tolerance TOL can be specified as an additional input argument.
%  
% See also PPLACE, STAB.

%	Author(s): M. Hromcik, M. Sebek 14-10-98
%	Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 1.0 $  $Date: 14-Oct-1998 10:28:34   $
%                         $Date: 04-Jul-2001  J.Jezek   $

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

eval('A = pol(A); B = pol(B);', 'error(peel(lasterr));');
[tv,var,A,B] = testvp(A,B);
if ~tv,
   warning('Inconsistent variables.');
end;

if strcmp(var,'s') | strcmp(var,'p'),
  error('A deadbeat controller only exists for discrete-time systems.');
end;  
      
switch opt,
  case 'l',
    if strcmp(var,'d') | strcmp(var,'z^-1'),
      if nargout <= 2,
	     eval('[Y,X] = pplace(B,A, pol(eye(A.s)), ''l'',tol);',...
	          'error(lasterr)');
	 
      else	
        eval('[Y,X,D,E,degT] = pplace(B,A, pol(eye(A.s)), ''l'',tol);',...
             'error(lasterr)');
        degT = [];
      end;	
    else,	% var = 'z' or 'q'
      if nargout <= 2,
        eval('[Y,X] = pplace(B,A, [0], ''l'',tol);',...
             'error(lasterr)');
      else 
        eval('[Y,X,D,E,degT] = pplace(B,A, [0], ''l'',tol);',...
             'error(lasterr)');
      end; 
    end;        	
	
  case 'r',
    if strcmp(var,'d') | strcmp(var,'z^-1'),
      if nargout <= 2,
         eval('[Y,X] = pplace(B,A, pol(eye(A.s)), ''r'',tol);',...
              'error(lasterr)');
      else	
      	eval('[Y,X,D,E,degT] = pplace(B,A, pol(eye(A.s)), ''r'',tol);',...
      	     'error(lasterr)');
      	degT = [];
      end;	
    else,	% var = 'z' or 'q'
      if nargout <= 2,
        eval('[Y,X] = pplace(B,A, [0], ''r'',tol);',...
         'error(lasterr)');
      else 
        eval('[Y,X,D,E,degT] = pplace(B,A, [0], ''r'',tol);',...
        'error(lasterr)');
      end; 
    end;        	

  otherwise 
    	error('Invalid command option.');
end;

if isempty(var), X.v = 'z'; Y.v = 'z'; D.v = 'z'; E.v = 'z'; end;

%end .. debe
