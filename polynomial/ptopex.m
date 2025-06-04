function [X,varargout] = ptopex(varargin)
%PTOPEX    Extreme polynomials for a polytope of polynomials
%
% The command
%      [AEX1,AEX2,...,AEXM] = PTOPEX(A0,A1,A2,...,AN,QBOUNDS)
% or                      X = PTOPEX(A0,A1,A2,...,AN,QBOUNDS), 
% create M = 2^N extreme polynomials AEX1,AEX2,...,AEXM 
% (or X = [AEX1;AEX2;...;AEXM]), for the polytope of polynomials
%      A = A0 + Q1*A1 + Q2*A2 + ...+ QN*AN
% defined by polynomials A0,A1,A2,...,AN and bounds Q1,Q2,...,QN.
% The uncertainty intervals [QiMIN,QiMAX] are given by
%      QBOUNDS = [ Q1MIN, Q1MAX | 
%                | Q2MIN, Q2MAX | 
%                |      ...     | 
%                | QNMIN, QNMAX ]
%
% The extreme polynomials are computed by enumerating A{Q} in the
% extremes of the parameter bounding set such as A{Q1EX,...,QNEX)
% where each QiEX is either QiMIN or QiMAX.
%
% See also PTOPLOT.

%       Author(s):  S. Pejchova, M. Sebek 08-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 15-Oct-1998 15:52:34   $
%       $Revision: 3.0 $  $Date: 28-Aug-2000  S. Pejchova   $

ni=nargin; no=nargout;  A=[];  Aii=[]; nr=max([0,ni-2]);

if ~ni, error('Not enough input arguments.');     
elseif (no>2^nr),
   error('Too many output arguments.');   
elseif ni<3,
   X = varargin{1}; return;
end;
qbounds=varargin{ni};
[rq,cq]=size(qbounds);
if (~isa(qbounds,'double')) | (ndims(qbounds) > 2) | length(qbounds)<1 |(rq~=2&cq~=2),
   error('The bounds array must be a nonempty double matrix.');
end;
if (rq==2)&(cq~=2), qbounds=(qbounds)';  rq=cq; cq=2; end;
if rq~=nr,
   error('The size of the bounds array must correspond to the number of polynomial inputs.');
end;
eval('A=pol(varargin{1});',...
     'error(''The input argument is not a polynomial matrix.'');');
[rA,cA]=size(A);  
Q_lim=zeros(2^nr,nr);

for ii=1:nr, 
  eval('Aii=pol(varargin{ii+1});',...
     'error(''The input argument is not a polynomial matrix.'');');
  if any([rA,cA]-size(Aii)),
     error('Matrices of inconsistent dimensions.');
  end;
  A=[A;Aii];
  Q_lim(:,ii)=repmat([ones(2^(ii-1),1);zeros(2^(ii-1),1)],2^(nr-ii),1);
end;
Q_min=Q_lim*(diag(qbounds(:,1)));   Q_min=[zeros(2^nr,1),Q_min];
Q_min=kron(Q_min,eye(rA));
Q_max=(~Q_lim)*(diag(qbounds(:,2))); Q_max=[ones(2^nr,1),Q_max];
Q_max=kron(Q_max,eye(rA));
X = Q_min*A+Q_max*A;
if no>1,
   X_out=X; X=X_out(1:rA,:);
   for ii=1:no-1,
      varargout{ii}=X_out(ii*rA+1:(ii+1)*rA,:);
   end;
end;

%end .. ptopex
