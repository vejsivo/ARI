function H = h2norm(N,D, arg1,arg2);
%H2NORM  H2-norm of a stable polynomial matrix fraction
%
% The commands
%    H2NORM(N,D)
%    H2NORM(N,D,'l')
% compute the H2-norm of the stable rational matrix G = D\N with 
% N and D left coprime. The command
%    H2NORM(N,D,'r') 
% computes the H2-norm of the stable rational matrix G = N/D with 
% N and D right coprime.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
% 
% See also POL/HINFNORM, POL/NORM.

%       Author(s): M. Hromcik, M. Sebek 22-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 27-Oct-1998 10:28:34   $
%       $Revision: 3.0 $  $Date: 23-Sep-1999 $
%       $Revision: 3.1 $  $Date: 23-Nov-1999 $
%		  $ Version 3 - 22-Nov-2000, Martin Hromcik $
%       $Revision $  $Date: 02-Jan-2002  J.Jezek  $
%                    $Date: 08-Jul-2002  J.Jezek  $

global PGLOBAL;
tol = PGLOBAL.ZEROING;
opt = 'l';

switch nargin,
  case 1,
    D = 1;
  case 2,
  case 3,
    if isnumeric(arg1),
      tol = arg1;
    elseif ischar(arg1),
      opt = arg1;
    else error('Invalid command option.');
    end;  
  case 4,
    if isnumeric(arg1),
      tol = arg1;
    elseif ischar(arg1),
      opt = arg1;
    else error('Invalid command option.');
    end;  
    
    if isnumeric(arg2),
      tol = arg2;
    elseif ischar(arg2),
      opt = arg2;
    else error('Invalid command option.');
    end;  
 
end;  	% switch ...

if ~all(size(tol)==1) | ~isreal(tol) | ...
      tol<0 | tol>1,
   error('Invalid tolerance.');
end;
 
eval('N = pol(N);', 'error(peel(lasterr));');
eval('D = pol(D);', 'error(peel(lasterr));');
Ns = N.s; Ds = D.s;
Nv = N.v; Dv = D.v;

% Variables:
if ~strcmp(Nv,Dv), 
  if isempty(Nv) | isempty(Dv),
    var = [Nv,Dv];
  else
    var = PGLOBAL.VARIABLE; 
    warning('Inconsistent variables.');
  end;    
else
  var = Nv;
end;    

if isempty(var), var = PGLOBAL.VARIABLE;
end;

[th,h,N,D] = testhp(N,D,var);
if ~th, warning('Inconsistent sampling periods.');
end;

% Sizes:
if Ds(1) ~= Ds(2), error('Denominator matrix is not square.');
end;
switch opt,
  case 'l',
     if Ds(1) ~= Ns(1),
        if Ds(1)==1,
           D = D*eye(Ns(1));
           Ds(1) = Ns(1); Ds(2) = Ns(1);
        else
           error('Matrices of inconsistent dimensions.');
        end;   
     end;
  case 'r',
     if Ds(2) ~= Ns(2),
        if Ds(2)==1,
           D = D*eye(Ns(2));
           Ds(1) = Ns(2); Ds(2) = Ns(2);
        else
           error('Matrices of inconsistent dimensions.');
        end;
     end;
  otherwise
    error('Invalid command option.');
end;    
if Ns(1)==0 | Ns(2)==0 | all(all(N==0)),
   H = 0; return;
end;

% Primeness:
switch opt,
  case 'l',
    if ~isprime(horzcat(N,D),'l'),
      error('Matrices are not left coprime.');
    end;
  case 'r',    
    if ~isprime(vertcat(N,D),'r'),
      error('Matrices are not right coprime.');
    end;
end;	%switch    

if any(strcmp(var,{'s','p'})), 
  % ..........
  % Continuous
  % ..........
  
  % Test if proper and stable:
  if ~isproper(N,D,opt,'strictly') 
    warning('Matrix fraction is not strictly proper. H2 norm is infinite.');
    H = Inf;
    return;	

  elseif ~isstable(D),
    warning('Denominator matrix is unstable. H2 norm is infinite.');
    H = Inf;
    return;	

  else	% strictly proper and strictly stable system
    switch opt,
      case 'l',
        %[N,D] = lmf2rmf(N,D);
        if rank(lcoef(D,'row')) < Ds(1),
          [D,aux,U] = rowred(D,tol);
          N = U*N;
        end;
        D = D.';
        N = N.';
      case 'r',
        if rank(lcoef(D,'col')) < Ds(1),
          [D,aux,U] = colred(D,tol);
          N = N*U;
        end;
    end;    % switch
    
    % Gramian G:
    
    % STEP 1:
    degcol = deg(D,'col')-1;
    degcol(degcol<=0) = 1;
    C = axxab(D,N'*N,degcol);
    if ~isproper(C,D,'r'), [aux,C] = rdiv(C,D); end;
    if isnan(C), error('The computation of the Gramian failed.'); end;

    % STEP 2:
    Dl = lcoef(D,'col');
    % Cl:
    Cs = C.s;
    Cc = C.c;
    Cl = zeros(Cs);
    c = deg(D, 'col');
    for i=1:size(D,2),
      %coli = Cc(:,i);
      %coli_c = coli.c;
	   coli_c = Cc(:,i,:);      
      %eval('Cl(:,i) = C{c(i)-1}(:,i);', 'Cl(:,i) = zeros(Cs(1),1);');
      eval('Cl(:,i) = coli_c(:,:,c(i));', 'Cl(:,i) = zeros(Cs(1),1);');
    end;

    G = Cl/Dl;
    if G(1,1)<0, G = -G; end;
    
    % H2 norm:
    H = sqrt(trace(G));
  end;  

else 
  % ........
  % Discrete
  % ........
  
  % Test if stable:
  if ~isstable(D) & D.d >= 1,
    warning('Denominator matrix is unstable. H2 norm is infinite.');
    H = Inf;
    return;
  end;  
  
  if any(strcmp(var,{'d','z^-1'})),
    [N,D] = reverse(N,D,opt,tol);
    N.v = 'z'; D.v = 'z';
    var = 'z';
  end;

  % Make proper:
  while ~isproper(N,D,opt),
    D = shift(D,1,var);
  end;  

  switch opt,
      case 'l',
        %[N,D] = lmf2rmf(N,D);
        if rank(lcoef(D,'row')) < Ds(1),
          [D,aux,U] = rowred(D,tol);
          N = U*N;
        end;
        D = D.';
        N = N.';
      case 'r',
        if rank(lcoef(D,'col')) < Ds(1),
          [D,aux,U] = colred(D,tol);
          N = N*U;
        end;
  end;    % switch

  P = axxab(D,N'*N); 
  if ~isproper(P,D,'r'), 
    [aux,P] = rdiv(P,D);
  end; 
      
  Dl = lcoef(D,'col');
  % Pl:
  Ps = P.s;
  Pc = P.c;
  Pl = zeros(Ps);
  c = deg(D, 'col');
  for i=1:size(D,2),
    %coli = P(:,i);
    %coli_c = coli.c;
    coli_c = Pc(:,i,:);
    %eval('Pl(:,i) = P{c(i)}(:,i);', 'Pl(:,i) = zeros(Ps(1),1);');
    eval('Pl(:,i) = coli_c(:,:,c(i)+1);', 'Pl(:,i) = zeros(Ps(1),1);');
  end; 
  gama0 = Pl/Dl;
  H = sqrt(trace(2*gama0));    
end;   

%end .. @pol/h2norm
