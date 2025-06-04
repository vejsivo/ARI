function y = isproper(N,D,varargin)
%ISPROPER  Test if polynomial matrix fraction is proper
%
% For polynomial matrices N,D, the commands
%    ISPROPER(N,D)
%    ISPROPER(N,D,'l') 
% return 1 if the left matrix fraction D\N is a proper rational 
% function, and 0 if it is not. The command
%    ISPROPER(N,D,'r') 
% returns 1 if the right matrix fraction N/D is a proper rational 
% function, and 0 otherwise. The command
%    ISPROPER(N,D,'strictly'), 
% possibly combined with the 'l' or 'r' option, tests if the matrix 
% fraction is strictly proper. The string 'strictly' may be shortened 
% to 'strict', 'str' or 's'.
%
% If D is omitted or empty then D = 1  is taken.
%
% For fraction in 's','p','z','q', the properness means the
% behaviour in point Infinity. For fraction in 'z^-1' or 'd',
% in point 0.
%
% A tolerance TOL may be specified as an additional (third, fourth
% or fifth) input argument. Its default value is the global zeroing
% tolerance.
%
% See also ISPRIME.

%       Author(s): M. Hromcik, M. Sebek 22-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 1.0 $  $Date: 22-Oct-1998 10:28:34   $
%       $          3.0 $  $Date: 02-Aug-2000  J.Jezek  $            
% 	     $          3.1 $  $20-Oct-2000, M. Hromcik - one argument accepted$ 
%		  $			 3.2 $  $22-Nov-2000,  bug fix (z^-1 & d with 's' option) $
%                         $02-Jul-2001, J. Jezek - logical output, error msgs $

global PGLOBAL;

tol = PGLOBAL.ZEROING; 
opt = 'l';
strict = 0;
twoarg = 0;

eval('N = pol(N);', 'error(peel(lasterr));');

if nargin>=2 & ~isempty(D),
   if ischar(D),
      if D(1)=='s',
         strict = 1;
      else
         opt = D;
      end;
   else
      eval('D = pol(D);', 'error(peel(lasterr));');
      twoarg = 1;
   end;
end;

lv = length(varargin);
for i = 1:lv,
   arg = varargin{i};
   if isnumeric(arg) & all(size(arg))==1 & ...
         isreal(arg) & arg>=0 & arg<=1,
      tol = arg;
   elseif ischar(arg),
      if arg(1)=='s',
         strict = 1;
      else
         opt = arg;
      end;
   else
      error(['Invalid ',nth(i+2),' argument.']);
   end;
end;

if ~twoarg, D = pol(1);
end;
[tv,var,N,D] = testvp(N,D);
if tv~=1,
   warning('Inconsistent variables.');
end;
[th,h,N,D] = testhp(N,D,var);
if th==0,
   warning('Inconsistent sampling periods.');
end;

switch opt,

  case 'l',	% D\N

    % Size consistency:
    eval('[N,D] = testdnd(N,D,opt);', ...
      'error(peel(lasterr));');

    % D must have full rank:
    if rank(D) < size(D,1),
      warning('Denominator matrix is singular.');
      y = logical(0);
      return;
    end;
    
    Nd = N.d;  % empty or zero
    if isempty(Nd) | Nd<0,
       y = logical(1); return;
    end;
    
    if any(strcmp(var, {'d','z^-1'})),
      td = min(tdeg(N),tdeg(D));
      Dc = D.c;
      y = ( rank(Dc(:,:,td+1)) == size(D,1) );
      if strict,
        Nc = N.c;
        y = y & all(all(Nc(:,:,td+1)==0));
      end;

    else	% the other vars    
      % Reduce D if necessary:
      Dlr = lcoef(D,'row');
      if rank(Dlr) < min(D.s),
        [D,aux,U] = rowred(D,tol);
        N = U*N;
      end;
    
      if strict,
        y = all( deg(N,'row') < deg(D,'row') );       
      else
        y = all( deg(N,'row') <= deg(D,'row') );
      end;
    end;
    
    
  case 'r',	% N/D

    % Size consistency:
    eval('[N,D] = testdnd(N,D,opt);', ...
      'error(peel(lasterr));');

    % D must have full rank:
    if rank(D) < size(D,1),
      warning('Denominator matrix is singular.');
      y = logical(0);
      return;
    end;
    
    Nd = N.d;  % empty or zero
    if isempty(Nd) | Nd<0,
       y = logical(1); return;
    end;
    
    if any(strcmp(var, {'d','z^-1'})),
      Dc = D.c; 
      ptr = 0;
      while all( Dc(:,:,ptr+1) == 0 ), ptr = ptr+1;D = shift(D,-1); end;
      
      N = shift(N, -ptr);
      Dc = D.c;
      y = ( rank(Dc(:,:,1)) == size(D,1) ) & ~isa(N,'tsp');
      if strict,
        N = tsp(N); 
        y = y & all(all(Nc(:,:,1)==0));
      end;

    else 	% the other vars   
    
      % Reduce D if necessary:
      Dlc = lcoef(D,'col');
      if rank(Dlc) < min(D.s),
        [D,aux,U] = colred(D,tol);
        N = N*U;
      end;
    
      if strict,
        y = all( deg(N,'col') < deg(D,'col') );       
      else
        y = all( deg(N,'col') <= deg(D,'col') );
      end;   
    end;
    
  otherwise error('Invalid command option.');

end;
    
%end .. @pol/isproper    
    