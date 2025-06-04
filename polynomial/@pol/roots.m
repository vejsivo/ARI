function r = roots(A,varargin);
%ROOTS  Roots of a polynomial matrix 
%
% The command
%    ROOTS(P) 
% computes the roots of a polynomial matrix P. If P is square
% nonsingular then the command
%    ROOTS(P,METHOD) 
% allows the user to specify the method that is used:
%    'det': computes the roots as the roots of the determinant of P.
%	  'eig': computes the roots as the generalized eigenvalues
%	         of the block companion matrix corresponding to P. 
%
% The commands
%    ROOTS(P,'all')
%    ROOTS(P,'det','all') 
% compute the infinite roots of P as well. If P is square then
%    ROOTS(P, 'eig', 'all') 
% computes the generalized eigenvalues of the block companion matrix 
% corresponding to P and leaves all large roots and all roots computed 
% as Inf and NaN in place.
%
% For polynomials in 'z' or in 'z^-1', it is possible to prescribe
% the variable for computing roots, the default being P.v . So, e.g.
%    ROOTS(z-0.5)  or  ROOTS(z-0.5,'z')  or  ROOTS(1-0.5*zi,'z')
% yields 0.5, whereas
%    ROOTS(1-0.5*zi)  or  ROOTS(1-0.5*zi,'zi')  or  ROOTS(z-0.5,'zi')
% yields 2.
%
% The commands
%    ROOTS(P,TOL)
%    ROOTS(P,METHOD,TOL)
%    ROOTS(P,METHOD,'all',TOL) 
% allow the user to define the tolerance TOL for rank testing. 
% Moreover, roots whose absolute value is greater than 1/TOL are 
% supposed to be infinite. The default value for TOL is PGLOBAL.ZEROING.
%
% See also DET, RANK.

%       Author(s): M. Hromcik, M. Sebek 20-5-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 5.0 $  $Date: 6-Oct-1998 14:28:34 - zeroing of determinant added for 'det' method (line 134) $
%       $Revision: 6.0 $  $Date: 10-Oct-2000 10:30:00 - Version 3.0, M. Hromcik $
%                         $Date: 17-Jul-2001 J. Jezek   arg checking  $
%                         $Date: 08-Jul-2002 J. Jezek   z,z^-1        $
%                         $Date: 21-Nov-2002 J. Jezek   z,z^-1        $
%                         $Date: 28-Feb-2003 J. Jezek   warning       $

global PGLOBAL;

eval('A = pol(A);','error(peel(lasterr));');
var = A.v;
allr = 0;
met = 'det';
tol = PGLOBAL.ZEROING;
recip = 0;

ni = nargin;
if ni>=2,
   for i = 2:ni,
      arg = varargin{i-1};
      if ~isempty(arg),
         if isa(arg,'pol'),
            [vs1,vs2,vd] = size(arg);
            if all([vs1,vs2,vd]==1) & all(arg.c(:,:)==[0 1]),
               arg = arg.v;
            else
               error(['Invalid ',nth(i),' argument.']);
            end;
         end;
         if ischar(arg),
            if strcmp(arg,'all'),
               allr = 1;
            else
               if strcmp(arg,'zi'), arg = 'z^-1';
               end;
               I = strmatch(arg,{'s';'p';'z^-1';'d';'z';'q'},'exact');
               if ~isempty(I),
                  if ~isempty(var),
                     if ~strcmp(var,arg),
                        if (strcmp(var,'z') & strcmp(arg,'z^-1')) | ...
                              (strcmp(var,'z^-1') & strcmp(arg,'z')),
                           recip = 1;
                        else
                           error('Invalid variable for roots.');
                        end;
                     end;
                  end;
               else
                  met = arg;
               end;
            end;
         elseif isnumeric(arg),
            tol = arg;
         else
            error(['Invalid ',nth(i),' argument.']);
         end;
      end;
   end;
end;

if ~length(tol)==1 | ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

Asc = A.s;
As = Asc(1);
Ac = A.c;
Ad = A.d;

if isempty(Ac) | Ad < 1,	% constant or zero or empty
  r = zeros(0,1);
  return;
end;

if recip,
   argallr = allr; allr = 1;
end;

rankA = rank(A);
is_sq_and_reg = ~any(Asc-As) & rankA == min(Asc);

if ~is_sq_and_reg			
% .............................
% Non square or singular matrix 
% .............................

  if ~strcmp(PGLOBAL.VERBOSE, 'no'), 
    warning(['Matrix is singular or not square.',...
             ' Possible METHOD option is not applied.']);
  end;  
  % Square down P (using rand) twice and choose genuine zeros: 
  d1 = det( rand(rankA, As)*A*rand(Asc(2), rankA), tol );
  %d1 = pzer(d1, tol*max(abs(d1.c(:))));		% !!!!!!!!!!!!
  r = roots(flipdim(d1.c,3));
  d2 = det( rand(rankA, As)*A*rand(Asc(2), rankA), tol );
  d2c = flipdim(d2.c,3);
  d2c = d2c(:).';
  for i = 1:length(r),
     if abs(polyval(d2c, r(i))) > tol, r(i) = NaN;
     end;
  end;   

else
% .......................
% Square full rank matrix
% .......................
 
 switch met,
  
  case 'det',
    detA = det(A,'fft', tol);

    if isinf(detA.d), 
      if ~strcmp(PGLOBAL.VERBOSE, 'no'), 
        warning(['Matrix has a zero determinant to working precision.', ...
                 'Finite roots computed via ''eig'' method.']);
      end; 
      r = roots_eig(Ad,As,Ac);
    else
      %detA = pzer(detA, tol*max(abs(detA.c(:))));	% !!!!!!!!!!!!!!!
      r = roots(flipdim(detA.c,3)); 
    end;  
  
  case 'eig',
    r = roots_eig(Ad,As,Ac);
  
  otherwise
    error('Invalid command option.');

 end;	% switch

end;

% ...........
% ALL option
% ...........
if ~allr,
  % Remove NaN, Inf and large roots
  r = r(~isnan(r));
  sw = warning; warning off;
  r = r(abs(r)<1/tol); warning(sw);
  if isempty(r), r = zeros(0,1);
  end;
  return;

else	% ALL
  if is_sq_and_reg & strcmp(met, 'eig'),	
    % EIG method and square regular A - leave the result unchanged
  else
    % Remove NaN, Inf and large roots and add infinites
    r = r(~isnan(r));
    sw = warning; warning off;
    r = r(abs(r)<1/tol); warning(sw);
    if isempty(r), r = zeros(0,1);
    end;
    [N,D] = reverse(A, eye(As), tol);
  
    % ............................................................
    % Compute the number of Inf roots of A as the zero roots of N:
  
    if ~is_sq_and_reg
      % A is singular and/or non square - random square-down:
      d3 = det( rand(rankA, As)*N*rand(Asc(2),rankA), tol );
      %d3 = pzer(d3, tol*max(abs(d3.c(:))));	% !!!!!!!!!!!!!!!!!!!!!!!!
      rinf = roots(flipdim(d3.c,3));
      rinf = rinf(abs(rinf) < tol);
      d4 = det( rand(rankA, As)*N*rand(Asc(2), rankA), tol );
      d4c = flipdim(d4.c,3);
      d4c = d4c(:).';
      for i = 1:length(rinf),
        if abs(polyval(d4c, rinf(i))) > tol, rinf(i) = NaN; end;
      end;   
      rinf = rinf(~isnan(rinf));
      rinf = Inf * ones(size(rinf));
      if isempty(rinf), rinf = zeros(0,1);
      end;
      r = [r; rinf];    

    else
      % Square non singular matrix A, MET ~= EIG:
      detN = det(N,'fft');
      if isinf(detN.d), 
        if ~strcmp(PGLOBAL.VERBOSE, 'no'), 
           warning(['Matrix has a zero determinant to working precision.', ...
                    'Finite roots computed via ''eig'' method.']);
        end; 
        rinf = roots_eig(N.d,N.s(1),N.c);
      else
        %detN = pzer(detN, tol*max(abs(detN.c(:))));	% !!!!!!!!!!!!!!!!!!!!!!!
        rinf = roots(flipdim(detN.c,3)); 
      end;  
      rinf = rinf(abs(rinf)<tol);	% Nearly zero roots
      rinf = Inf * ones(size(rinf));
      r = [r; rinf];
    end;  
  end;  
end; 	%if ~all ...  

if recip,
   if ~argallr,
      if ~isempty(r),
         zeroi = find(r==0);
         r(zeroi) = [];
      end;
   end;
   saved_w = warning; warning off;
   r = 1./r;
   warning(saved_w);
end;

% .......................................................................
% sub-function for 'eig' roots computation of square non singular matrix: 
function r = roots_eig(Ad, As, Ac);
  global PGLOBAL;
  C = diag(ones(1,As*(Ad-1)), As);
  Acrow = Ac(:,:);
  C( (Ad-1)*As+1:Ad*As, :) = -Acrow(:,1:As*Ad);
  D = eye(As*Ad);
  D(As*(Ad-1)+1:As*Ad,As*(Ad-1)+1:As*Ad) = Ac(:,:,Ad+1);
  r = eig(C,D);
%end .. roots_eig
% .......................................................................

%end .. @pol/roots
