function [Ak,Bk,Ck,Dk,Ek,gopt,clpoles] = dssrch(A,B,C,D,E,nmeas,ncon,gmin,gmax,p10,p11)
%DSSRCH  Search for the optimal solution of a descriptor H-infinity problem
%
% The command
%    [Ak,Bk,Ck,Dk,Ek,gopt,clpoles] = DSSRCH(A,B,C,D,E,nmeas,ncon,gmin,gmax[,pars][,'show'])
% computes the central H-infinity optimal solution for the standard plant with
% transfer matrix 
%   G(s) = C*(s*E-A)^{-1}*B + D, 
% ncon control inputs and nmeas measured outputs. 
% The optimal compensator has the transfer matrix 
%    K(s) = Ck*(s*Ek-Ak)^{-1}*Bk + Dk. 
% The output parameter gopt is the minimal H-infinity norm and clpoles 
% represents the closed-loop poles. 
% The input parameter gmin is the lower limit of the initial search interval
% and gmax the upper limit.
%
% Features:
% (1) If no solution exists within the initial search interval [gmin,gmax] then
%     the interval is automatically expanded.
% (2) The macro also computes optimal solutions for generalized plants that
%     have unstable fixed poles.
% (3) The macro computes both type 1 optimal solutions (gopt coincides with
%     the lowest value of gamma for which spectral factorization is possible)
%     and type 2 optimal solutions (gopt is larger than this lowest value.
% (4) When the search is completed the compensator is reduced to minimal
%     dimension. If a proper compensator exists then Ek is the unit matrix.
% (5) Initially the search is binary. When the solution is of type 2 and close
%     to optimal then a secant method is used to obtain fast convergence to
%     an accurate solution.
% (6) Before the search is initiated the descriptor representation of the
%     generalized plant is ``regularized'' in the sense that it is modified
%     so that D12 and D21 have full rank.
% (7) The optional input parameter pars is a vector with components
%        pars = [acc tolrho thrrho maxpole tolrnk tolstable]
%     with
%        acc  Accuracy. If gmax-gmin < acc then the search is terminated. The 
%             default value is 1e-4.
%        tolrho Tolerance on the (nonnegative) singularity parameter rho. A type 2 
%             optimal solution is achieved exactly if rho = 0. If rho < tolrho then
%             the search is terminated. The default value is 1e-8.
%        thrrho Threshhold on rho. If rho < thrrho then the search switches from
%             binary to secant. The default value is 0.01.
%        maxpole Largest size of an open- or closed-loop pole. If the size of a pole 
%             is greater than maxpole then the pole is considered to be a pole
%             at infinity. A fixed open-loop pole is considered unstable if its
%             real part is greater than -1/maxpole. The default value of maxpole 
%             is 1e+6.  
%        tolrnk Tolerance in the rank test to determine whether an open-loop pole
%             is uncontrollable or observable and, hence, a fixed pole. The default
%             value is 1e-8
%        tolstable Tolerance used to test stability. The default value is 1e-8.
%
% If pars has p entries, with p < 6, then these entries are assigned to the first p 
% parameters listed.
%
% If the optional input parameter 'show' is present then the progress of the search
% is displayed.
%
% See also DSSHINF.

%    Author H. Kwakernaak, 1998
%    Copyright 1998 by PolyX Ltd.
%    Modified by J. Jezek, Aug 2001, arg checking

% Tolerances
% pars = [acc tolrho thrrho maxpole]


% Checks

if nargin < 9
   error('Not enough input arguments.');
end;

show = 0;
pars = [];
if nargin == 10 | nargin == 11 
   if isstr(p10)
      if strcmp(p10,'show')
         show = 1;
      else
         error('Invalid command option.')
      end
   elseif isnumeric(p10) & isreal(p10)
      pars = p10;
   else
      error('Invalid 10th argument.');
   end
end
if nargin == 11 
   if isstr(p11)
      if strcmp(p11,'show')
         show = 1;
      else
         error('Invalid command option.')
      end
   elseif isnumeric(p11) & isreal(p11)
      pars = p11;
   else
      error('Invalid 11th argument.');
   end
end

acc = 1e-4; tolrho = 1e-8; thrrho = 0.01; maxpole = 1e6; 
tolrnk = 1e-8; tolstable = 1e-8;
if ~isempty(pars)
   switch length(pars)
      case 1, acc = pars(1);
      case 2, acc = pars(1); tolrho = pars(2);
      case 3, acc = pars(1); tolrho = pars(2); thrrho = pars(3);
      case 4, acc = pars(1); tolrho = pars(2); thrrho = pars(3); 
              maxpole = pars(4);
      case 5, acc = pars(1); tolrho = pars(2); thrrho = pars(3); 
              maxpole = pars(4); tolrnk = pars(5);
      case 6, acc = pars(1); tolrho = pars(2); thrrho = pars(3); 
              maxpole = pars(4); tolrnk = pars(5); tolstable = pars(6);
      otherwise, error('Invalid 10th or 11th argument; too long.');
   end
end
if acc<0 | acc>1 | tolrho<0 | tolrho>1 |tolrnk<0 | tolrnk>1 | ...
      tolstable<0 | tolstable>1
   error('Invalid tolerance.');
end

if ~isnumeric(A) | ndims(A)>2 | ~isnumeric(B) | ndims(B)>2 | ...
   ~isnumeric(C) | ndims(C)>2 | ~isnumeric(D) | ndims(D)>2 | ...
   ~isnumeric(E) | ndims(E)>2
   error('Invalid 1st - 5th argument; must be numerical matrices.');
end;
[nE,mE] = size(E); [nA,mA] = size(A); [nB,mB] = size(B); 
[nC,mC] = size(C); [nD,mD] = size(D);
if nE~=mE | nA~=mA | nE~=nA | mC~=nA | mA~=nB | nD~=nC | mD~=mB
   error('Matrices of inconsistent dimensions.')
end

% Initial checks

if ~isnumeric(nmeas) | length(nmeas)~=1 | ~isreal(nmeas) | ...
   ~isnumeric(ncon)  | length(ncon)~=1  | ~isreal(ncon)
   error('Invalid 6th or 7th argument.');
end;
if ~isnumeric(gmin) | length(gmin)~=1 | ~isreal(gmin) | ...
   ~isnumeric(gmax) | length(gmin)~=1 | ~isreal(gmin)   
   error('Invalid 8th or 9th argument.');
end;
if gmax < gmin
   error('gmax is less than gmin.')
end

% Initializations

global PGLOBAL;
eval('PGLOBAL.FORMAT;', 'painit;')
verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

rhomin = Inf; rhomax = Inf;
Gmin = gmin; Gmax = gmax;

% Regularize the data

eval('[A,B,C,Dnew,E] = dssreg(A,B,C,D,E,nmeas,ncon);', ...
   'error(peel(lasterr));');
if ~isempty(D) & any(any(Dnew~=D))
   if show
      fprintf('\ndssrch: plant transformation to regularize it\n')
   end
   D = Dnew;
end
k2 = ncon; m2 = nmeas; 
k1 = mB-k2; m1 = nC-m2; 
p = length(E);
B2 = B(:,k1+1:k1+k2);
C2 = C(m1+1:m1+m2,:);
D22 = D(m1+1:m1+m2,k1+1:k1+k2);

% Compute the fixed poles

olpoles = eig(A,E); 
olpoles = olpoles(find(~isnan(abs(olpoles))));
olpoles = olpoles(find(abs(olpoles)<maxpole));
fxpoles = [];
if ~isempty(olpoles)
   for i = 1:length(olpoles)
      s = olpoles(i);
      if rank([s*E-A B2],tolrnk) < p | rank([s*E-A; C2],tolrnk) < p
         fxpoles = [fxpoles; s];
      end
   end
end
unstablefxpoles = 0;
if ~isempty(fxpoles)
   fprintf('\nFixed poles\n\n')
   for i = 1:length(fxpoles)
      fprintf('     %g',fxpoles(i))
      if real(fxpoles(i)) >= -1/maxpole
         unstablefxpoles = unstablefxpoles + 1;
         fprintf('\tunstable\n')
      else
         fprintf('\n')
      end
   end
end
   

% Search for the optimal solution

if show
   fprintf('\nSearch for optimal solution\n')
   fprintf('\n     action\tgamma\n')
   fprintf('     ------\t-----\n')
end

status = 'testgmax';
notdone = 1;

while notdone
   Status = status;
   if strcmp(status,'testgmax')
      if gmax == Gmax*2^5
         fprintf('dssrch: gmax has been doubled 4 times without finding\n')
         fprintf('a stabilizing compensator at gamma = gmax = %g\n',gmax/2)
         error(' ')
      else
         gamma = gmax;
      end
   elseif strcmp(status,'testgmin')
      if gmin == Gmin/2^5
         fprintf('dssrch: gmin has been halved 4 times but a stabilizing\n')
         fprintf('compensator still exists at gamma = gmin = %g\n',2*gmin)
         error(' ')
      else
         gamma = gmax;
      end
      gamma = gmin;
   elseif strcmp(status,'binsearch')
      gamma = (gmin+gmax)/2;
   elseif strcmp(status,'secsearch')
       gamma = (gmax*rhomin+gmin*rhomax)/(rhomax+rhomin);
   end
   [Ak,Bk,Ck,Dk,Ek,clpoles,rho] = dsshinf(A,B,C,D,E,nmeas,ncon,gamma,'notest');
   if ~isempty(rho)
      stable = (length(find(real(clpoles)>=0)) == unstablefxpoles);
   end
   if strcmp(status,'testgmax')
      if isempty(rho) 
         gmax = 2*gmax;
      elseif ~stable
         gmax = 2*gmax;
      else
         rhomax = rho;
         status = 'testgmin';
      end
   elseif strcmp(status,'testgmin')
      if isempty(rho)
         status = 'binsearch';
      elseif stable
         gmin = gmin/2;
      else
         rhomin = rho;
         status = 'binsearch';
      end
   elseif strcmp(status,'binsearch') | strcmp(status,'secsearch')
      if isempty(rho)
         gmin = gamma; rhomin = rho;
      elseif stable
         gmax = gamma; rhomax = rho;
      else
         gmin = gamma; rhomin = rho;
      end
      if ~isempty(rho)
         notdone = ( (Gap > acc) & (rho > tolrho) );
      end
      if ~isempty(rhomin) & ~isempty(rhomax) & isfinite(rhomin) & isfinite(rhomax)
         if rho < thrrho
            status = 'secsearch';
         end
      end
   end     
   Gap = abs(gmax-gmin);
   if show
      fprintf('     %s\t%g\n',Status,gamma)
   end
end % notdone


% Reduce the compensator and complete the output

gopt = gmax;
[ak,bk,ck,dk,ek] = dsshinf(A,B,C,D,E,nmeas,ncon,gopt,'notest',1/maxpole);
[Ak,Bk,Ck,Dk,Ek] = dssmin(ak,bk,ck,dk,ek,1/maxpole);
Ec = [            E                zeros(size(E,1),size(Ek,2))
       zeros(size(Ek,1),size(E,2))               Ek           ];
Ac = [    A       B2*Ck
        Bk*C2 Ak+Bk*D22*Ck ];
clpoles = roots(pol([-Ac Ec],1));
     
%end .. dssrch
     