function Ps = prod(P, dim, tol);
%PROD   Element-wise product of polynomial
%
% For  vectors, PROD(X) is the product of the elements of X. 
% For matrices, PROD(X) is a row vector with the product over each column. 
%
% PROD(X,DIM) works along the dimension DIM.
%
% PROD(X,DIM,TOL) or PROD(X,[],TOL) works with zeroing specified by the
% input tolerance TOL .

%       Author(s): M. Hromcik, M. Sebek 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date:  7-May-1998 13:28:34   $
%       $Revision: 3.0 $  $Date: 13-Oct-1999 12:00:00   J. Jezek  $
%                         $Date: 03-Feb-2000 12:00:00   J. Jezek  $
%                         $Date: 26-May-2000 15:00:00   J. Jezek  $
%                         $Date: 28-Feb-2003            J. Jezek  $

global PGLOBAL;

eval('P = pol(P);','error(peel(lasterr))');

ni = nargin;
if ni==1,
  dim = []; tol = PGLOBAL.ZEROING;
elseif ni==2 | isempty(tol),
  tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
     error('Invalid tolerance.');
  end;
end;

Pc = P.c; 
Psi = P.s; 
Pd = P.d;
Pv = P.v;
Ph = P.h;

if isempty(dim),
   if isempty(Pd),
      Ps = pol(prod(zeros(Psi)));
      return;
   end;
   dim = find(Psi~=1);
   if isempty(dim),
      dim = 1;
   else
      dim = dim(1);
   end;
end;
if ~isa(dim,'double') | length(dim)~=1  | (dim~=1 & dim~=2),
   error('Invalid dimension.');
end;

if isempty(Pd),
   Ps = pol(prod(zeros(Psi),dim));
   return;
end;

switch dim,
 case 1,
  Ps = ones(1,Psi(2));
  Pi.d = Pd; Pi.s = [1 Psi(2)]; 
  for i = 1:Psi(1),
    Pi.c = Pc(i,:,:); Pi.v = Pv;
    Pi.h = Ph; Pi.u = []; Pi.version = 3.0;
    Pip = class(Pi, 'pol');
    Ps = times(Ps, Pip, 0);
  end;

 case 2
  Ps=ones(Psi(1), 1);
  Pi.d = Pd; Pi.s = [Psi(1) 1]; 
  for i=1:Psi(2),
    Pi.c = Pc(:,i,:); Pi.v = Pv;
    Pi.h = Ph; Pi.u = []; Pi.version = 3.0;
    Pip = class(Pi, 'pol');
    Ps = times(Ps, Pip, 0);
 end;
 
end;	% switch

% zeroing:
me = min(abs(nonzeros(Pc)));
Ps.v = Pv;  Ps.h = Ph;
if ~isempty(me),
  Ps = pzer(Ps,tol*me);
end;

%end .. @pol/prod
	
  