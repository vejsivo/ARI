function [Ne,De] = reverse(N,D, arg3, arg4);
%REVERSE   Reverse the variable of a polynomial matrix 
%           or a polynomial matrix fraction
%
% For a single polynomial matrix P = P0 + P1*v + ... + Pn*v^n, 
% the command 
%     PR = REVERSE(P)
% reverses the order of particular coefficients of P, namely 
% PR = P0*v^n + P1*v^(n-1) + ... + Pn.
%
% Given the polynomial matrices P and Q the commands
%    [N,D] = REVERSE(P,Q)
%    [N,D] = REVERSE(P,Q,'l') 
% compute the polynomial matrices N and D such that 
%    INV(D(VAR)) * N(VAR) = INV(Q(1/VAR)) * P(1/VAR)
% If P and Q are left coprime then so are N and D. VAR is the variable
% of P and Q. The command
%    [N,D] = REVERSE(P,Q,'r') 
% computes the matrices N and D such that 
%    N(VAR) * INV(D(VAR)) = P(1/VAR) * INV(Q(1/VAR))
% If P and Q are right coprime then N and D are right coprime as well.
%
% A tolerance TOL may be specified as an additional input argument.
%
% If 
%    T(VAR) = INV(Q(VAR)) * P(VAR) 
% is the transfer matrix of a linear discrete-time system in the forward 
% shift operator (usually denoted as 'z') then
%    [N, D] = REVERSE(P, Q)  
% computes the left coprime numerator and denominator matrices of the 
% transfer matrix T expressed in the backward shift operator (usually 
% denoted as '1/z', 'z^-1' or 'd'). The inverse transformation may be 
% performed similarly.
%
% See also LMF2SS, RMF2SS, LMF2RMF.

%	Author(s): M. Hromcik, M. Sebek 20-9-98
%	Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 28-Sep-1998 10:28:34          $
%       $Revision: 3.0 $  $Date: 3-10-2000 - one-input added   $
%                         $Date; 19-Jul-2001 J.Jezek sampl per $
%                         $Date: 21-Nov-2002 J.Jezek bug       $ 

% 'l' .. D^-1 * N; rowred([D N]);
% 'r' .. N * D^-1; colred([D;N]);

global PGLOBAL;

eval('N = pol(N);', 'error(peel(lasterr));');

tol = PGLOBAL.ZEROING;
opt = 'l';

if nargin == 1,
   %error('Not enough input arguments.');
   Ne = N;
   eval('Ne.c = flipdim(N.c,3);', '');
   Ne = pzer(Ne,tol);    %  J.Jezek  21-Nov-2002 
   if nargout>1,
      error('Too many output arguments.');
   end;
   return;
end;

eval('D = pol(D);', 'error(peel(lasterr);');

if isempty(N.c) & isinf(N.d), N.c = zeros(N.s); N.d = 0; end;
if isempty(D.c) & isinf(D.d), D.c = zeros(D.s); D.d = 0; end;

switch nargin,
  
  case 1,
    
  case 2,
  
  case 3,
     if isnumeric(arg3),
        if ~isempty(arg3), tol = arg3;
        end;
    elseif isa(arg3,'char'), opt = arg3;
    else error('Invalid command option.');
    end;
  
  case 4,
    if isnumeric(arg3),
        if ~isempty(arg3),tol = arg3;
        end;
    elseif isa(arg3,'char'), opt = arg3;
    else error('Invalid command option.');  
    end;
    if isnumeric(arg4),
       if ~isempty(arg4), tol = arg4;
       end;
    elseif isa(arg4,'char'), opt = arg4;
    else error('Invalid command option.');
    end;    

end;    % switch

if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

Nsc = N.s;
Ns1 = Nsc(1);
Ns2 = Nsc(2);
Dsc = D.s;
Ds = Dsc(1);
if Ds ~= Dsc(2), error('Denominator matrix is not square.'); end;
if rank(D)<Ds,
 warning('Denominator matrix is singular to working precision.'); 
end;

switch opt
  case 'l',
    [Ne, De] = lmf_dz(N,D,Ds,Dsc,Ns1, Nsc, tol);
  case 'r',
    [Ne, De] = lmf_dz(N.',D.',Ds,Dsc,Ns2, fliplr(Nsc), tol); 
    Ne = Ne.';
    De = De.';
  otherwise error('Invalid command option.');
end;  

% Sub function for LMF:
function [Ne, De] = lmf_dz(N, D, Ds, Dsc, Ns1, Nsc, tol);

global PGLOBAL;

if Ds ~= Ns1, 
     if all(Dsc == 1),
       D = eye(Ns1)*D;
       Dsc = D.s;
       Ds = Dsc(1);
     else,  
       error('Matrices of inconsistent dimensions.'); 
     end;
end;

[tv,var,N,D] = testvp(N,D);
if tv~=1,
   warning('Inconsistent variables.');
end;
  
[th,h,N,D] = testhp(N,D,var);
if ~th,
   warning('Inconsistent sampling periods.');
end;

DN = [D N];
degrowDN = deg(DN, 'row');
if rank(lcoef(DN,'row')) < Ds, 
  DN = rowred(DN, tol);
  degrowDN = deg(DN, 'row');
end;
degrowDN(isinf(degrowDN)) = 0;
   
DNs = DN.s;
DNc = DN.c;
DNd = DN.d;
Aux = zeros(size(DNc));

for i = 1:DNs(1),
   	Auxi = DNc(i,:,1:degrowDN(i)+1);
  	Aux(i,:,1:degrowDN(i)+1) = flipdim(Auxi, 3);
end;


%De.d = DNd;
%De.s = Dsc;
%De.c = Aux(:,1:Ds,:);
%De.v = var;
%De.u = [];
%De.version = 2.0;
%De = class(De,'pol');
Dec = Aux(:,1:Ds,:);
De = pol(Dec(:,:),DNd,var);
De.h = h;

%Ne.d = DNd;
%Ne.s = Nsc;
%Ne.c = Aux(:,Ds+1:end,:);
%Ne.v = var;
%Ne.u = [];
%Ne.version = 2.0;
%Ne = class(Ne,'pol');
Nec = Aux(:,Ds+1:end,:);
Ne = pol(Nec(:,:),DNd,var);
Ne.h = h;

me = min(abs([nonzeros(N.c); nonzeros(D.c)]));
if ~isempty(me),
  Ne = pzer(Ne, tol*me);
  De = pzer(De, tol*me);
end;  
 
%end .. lmf_dz

%end .. @pol/reverse
