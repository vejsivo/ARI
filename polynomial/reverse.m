function [Ne,De] = reverse(N,D, arg3, arg4)
%REVERSE    Reverse the variable of a constant matrix
%             or a constant matrix fraction
%
% For a single constant matrix P, the command
%     PR = REVERSE(P)    returns  PR = P .
%
% Given the constant matrices P and Q, the commands
%    [N,D] = REVERSE(P,Q)
%    [N,D] = REVERSE(P,Q,'l')
% return  N,D  such that  INV(D)*N = INV(Q)*P .
%
% The command
%    [N,D] = REVERSE(P,Q,'r')
% returns  N,D  such that  N*INV(D) = P*INV(Q) .
%
% This maco exists only for completeness. For more information
% and for possible further arguments, see POL/REVERSE.

%    Author: J.Jezek, 13-Aug-2001
%    Copyright(c) 2001 by Polyx, Ltd.

ni = nargin;
if ni < 1,
   error('Not enough input arguments.');
end;
eval('N = pol(N);', 'error(peel(lasterr));');
if ni == 1,
   if nargout == 1,
      Ne = reverse(N);
   else
      error('To many output arguments.');
   end;
else
   eval('D = pol(D);', 'error(peel(lasterr));');
   if ni == 2,
      eval('[Ne,De] = reverse(N,D);', ...
         'error(peel(lasterr));');
   else
      if ~(isa(arg3,'double') | isa(arg3,'char')),
         error('Invalid 3rd argument.');
      end;
      if ni == 3,
         eval('[Ne,De] = reverse(N,D,arg3);',...
            'error(peel(lasterr));');
      else
         if ~(isa(arg4,'double') | isa(arg4,'char')),
            error('Invalid 4th argument.');
         end;
         eval('[Ne,De] = reverse(N,D,arg3,arg4);', ...
            'error(peel(lasterr));');
      end;
   end;
end;

%end .. reverse

         