function [Q,R,Z] = hqr(A,J,tol)
% Given a square matrix A and a diagonal signature matrix J, the command
%
%   [Q,R,Z] = HQR(A,J)
%
% computes a hyperbolic QR factorization of A such that matrix
%
%   Q*A = R
%
% has as many zero rows as the rank defect of A,
% and Q a J-orthogonal matrix, i.e.
%
%   Q'*J*Q = J
%
% Vector Z contains indices of non-zero rows in matrix R,
% hence length(Z) = rank(A)
%
% The reduction is possible if and only if rank(A) = rank(A'*J*A).
% Otherwise, the function returns empty matrices Q and R
%
% An optional tolerance parameter can be specified. It is used as an
% absolute threshold to detect zero elements during the reduction process.
% By default it is the global zeroing tolerance
%
% See also: QR
  
% The algorithm consists in zeroing iteratively entries in each column
% of A by application of orthogonal (plane) rotations, row permutations
% and J-orthogonal (hyperbolic) rotations, see [D. Henrion, P. Hippe.
% Hyperbolic QR factorization for J-spectral factorization of
% polynomial matrices. LAAS-CNRS Research Report, Toulouse, France,
% February 2003]

% Written by D. Henrion and P. Hippe, November 22, 2002.
% Last modified by D. Henrion, February 10, 2003.
% Copyright 2002-2003 by PolyX, Ltd

if nargin < 3
 global PGLOBAL;
 if exist(PGLOBAL.ZEROING)
  tol = PGLOBAL.ZEROING;
 else
  tol = 1e-8;
 end
end;

[n,m] = size(A);
if m ~= n
  error('Input matrix must be square')
end
R = A; Q = eye(n);
J = sign(diag(J)); % keep only diagonal entries
if any(J==0)
  error('J must be a signature matrix')
end

U = eye(n); % column operations
perm = 1:n; % row permutations

r = 0; c = 0;
while c < n
 c = c+1; % column index
 if max(abs(R(r+1:n,c))) >= tol % entries to be cancelled ?
  r = r+1; fail = (r < n); again = 0;
  while fail
   fail = 0; cancel = 0;
   r1 = r; r2 = r1;
   while r2 < n % cancel entries in the column
    r2 = r2+1; a1 = R(r1,c); a2 = R(r2,c);
    if abs(a2) >= tol % non-zero a2
     if J(r1) == J(r2) % same sign, apply orthogonal rotation
      cancel = 1;
      t = abs(a1)+abs(a2);
      if abs(t) > tol
       x = a1/t; y = a2/t; z = sqrt(x^2+y^2); x = x/z; y = y/z;
       R([r1 r2],:) = [x y;-y x]*R([r1 r2],:);
       Q([r1 r2],:) = [x y;-y x]*Q([r1 r2],:);
      end
     else % different sign, apply hyperbolic rotation..
      da = abs(a1)-abs(a2);
      if abs(da) >= tol % ..when possible
       cancel = 1;
       if da < 0 % permute rows if necessary
        ind = 1:n; ind(r1) = r2; ind(r2) = r1;
        Q = Q(ind,ind); R = R(ind,:); J = J(ind);
        perm = perm(ind); a3 = a1; a1 = a2; a2 = a3;
       end
       t = a2/a1; x = 1/sqrt(1-t^2); y = x*t;
       R(r1,:) = x*R(r1,:)-y*R(r2,:); R(r2,:) = -y/x*R(r1,:)+R(r2,:)/x;
       Q(r1,:) = x*Q(r1,:)-y*Q(r2,:); Q(r2,:) = -y/x*Q(r1,:)+Q(r2,:)/x;
      else % no hyperbolic rotation
       fail = 1;
       if again % try to combine columns..
        d = c;
        while fail & (d < n)
         d = d+1; b1 = R(r1,d); b2 = R(r2,d);
         if abs(a1*b2-a2*b1) >= tol % .. when possible
          x = rand; y = sqrt(1-x^2);
          R(:,[c d]) = R(:,[c d])*[x y;-y x];
   	  U(:,[c d]) = U(:,[c d])*[x y;-y x];
    	  fail = 0;
         end;
        end
	if ~fail % and try same column again
	 r2 = r1; 
	end
       end
      end
     end
    end     
   end
   if fail % some hyperbolic rotation failed
    if cancel % but some entries were cancelled
     again = 1; % so try again
    else % and no entries were cancelled
     if again % already tried once so give it up
      Q = []; R = []; Z = []; return
     else
      again = 1; % try again
     end
    end
   end
  end
 end
end
% post-processing
Q(perm,perm) = Q; R(perm,:) = R*U'; Z = perm(1:r);



