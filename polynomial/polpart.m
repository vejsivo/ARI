function [Bl, Br, offset] = polpart(B, degree)
%POLPART  Polynomial matrix part extraction
%
% For square para-Hermitian polynomial matrix B, the command
%    BR = POLPART(B) 
% extracts the half part of B such that
%    B = BR' + BR                in continuous-time, and
%    B = BR' + SHIFT(BR,DEG(BR)) in discrete-time
%
% For square discrete-time polynomial matrix B, the command
%    [BL BR OFFSET] = POLPART(B,CONSTANT) 
% extracts left and right parts of B such that
%    B{CONSTANT} = BL{0} + BR{0}
%    B = BL' + SHIFT(BR,OFFSET)
% If CONSTANT is not specified then CONSTANT = FLOOR(DEG(B)/2).
%
% See also POL/CTRANSP.

%     Author: D. Henrion, January 19, 1999.
%     Copryright 1998-99 by Polyx, Ltd.
%     Modified by J. Jezek, Aug 2001, arg check

if nargin < 1,
   error('Not enough input arguments.');
end;
eval('B = pol(B);', 'error(peel(lasterr));');

if nargin == 1 | isempty(degree),
   degree = floor(deg(B) / 2);
else
   if ~isnumeric(degree) | length(degree)~=1 | ...
         ~isreal(degree) | degree<0,
      error('Invalid degree offset.');
   end;
end;

[n, cB] = size(B); degB = deg(B);

if n ~= cB,
 error('Matrix is not square.');
end;

typeB = B.var;
if strcmp(typeB, 's') | strcmp(typeB, 'p'),
 type = 'continuous';
else
 type = 'discrete';
end;

if (nargout == 2),
 if strcmp(type, 'continuous'),
  error('Too many output arguments for continuous-time polynomial matrix.');
 end;
end;

if (nargout < 2)
 residue = norm(B-B');
 if residue > eps,
  warning(['POLPART: Polynomial matrix is not para-Hermitian. Residue = ' num2str(residue)]);
 end;
end;

if strcmp(type, 'continuous'),

 Bl = zeros(n, n*(degB+1));
 for i = 0:2:degB,
  Bl(:,1+i*n:(i+1)*n) = B{i} / 2; % symmetric part
 end;

 for i = 1:2:degB,
  Bl(:,1+i*n:(i+1)*n) = tril(B{i}, -1); % anti-symmetric part
 end;

 Bl = pol(Bl, degB); pprop(Bl, typeB);

else % discrete-time

 if degree > degB,

  offset = degree; i = 0;
  while (i < degB) & abs(B{i} < eps),
   offset = offset - 1; i = i + 1;
  end;
  Bl = B';
  if nargout > 1,
   Br = pol(zeros(n));
  end;

 else  
 
  B{degree} = B{degree} / 2;
  Bl = pol(B{0:degree}, degree); pprop(Bl, typeB);

  offset = degree; i = 0;
  while (i < degree) & abs(Bl{i} < eps),
   offset = offset - 1; i = i + 1;
  end;
  Bl = Bl';

  if nargout > 1,
   Br = pol(B{degree:degB}, degB-degree); pprop(Br, typeB);
  end;

 end;

end;

%end .. polpart

