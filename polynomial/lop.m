function P = lop(varargin)
%LOP Specify a polynomial matrix or convert a constant matrix to a polynomial matrix
%
% The commmand
%    P =  LOP(A,D) 
% creates a polynomial matrix object
%        P(v) = P0 + P1*v + P2*v^2 + .. + PD*v^d
% from the constant matrix A = [Pd,...,P2,P1,P0] containing 
% the matrix coefficients ordered according to descending powers of v.
% If A is a 3-dim constant Matlab matrix, and [sz1,sz2,sz3] = size(A), 
% P = LOP(A) returns the corresponding polynomial matrix of degree D. 
% If (sz3-1) < D, then
% and P0 = zeros(sz1,sz2), ..., P(D-1) = A(:,:,2), PD = A(:,:,1). 
% If (sz3-1) >= D  and  AA = A(:,:) then
% P0 = AA(:,(end-((sz2*sz3)/(D+1))+1):end),..., PD = AA(:,1:(sz2*sz3)/(D+1)).
%
% The commmand
%    P =  LOP(A) 
% returns P = A if A is already a polynomial matrix object. If A is a 
% (standard MATLAB) constant matrix, then it returns the corresponding 
% zero degree polynomial matrix object
%
% If A is a 3-dim constant Matlab matrix, and [sz1,sz2,sz3] = size(A), 
% P = LOP(A) returns the corresponding polynomial matrix of degree D = sz3-1.
% and  P0 = A(:,:,end),..., PD = A(:,:,1). 
%
% All the above syntaxes may be followed by PropertyValue arguments.
% (Type "help pprop" for details on assignable properties)
%   
% If A is zero and/or d = -Inf then P is a zero polynomial matrix.
% If A and/or d is empty then an empty polynomial matrix object is
% created.
%
% If not specified by a PropertyValue argument then the variable string 
% 'v' copies the current value of the global variable symbol. The string 
% can be
%  * the continuous-time operator: 's' (default) or 'p'
%  * the discrete-time forward shift operator: 'z' or 'q'
%  * the discrete-time backward shift (delay) operator: 'd' or 'z^-1'
%
% If A is a two-sided polynomial matrix then the commands
%    P = LOP(A),  P = LOP(A,'z'),  P = LOP(A,'z^-1')
% return the corresponding polynomial matrix whenever possible.
%
% If A is a symbolic matrix of the Symbolic Toolbox then
% P = LOP(A) converts A to a polynomial matrix P. The variable string
% copies the current value of the global variable symbol.
%
% If A is a 3-dim constant Matlab matrix, and [sz1,sz2,sz3] = size(A), 
% P = LOP(A) returns the corresponding polynomial matrix of degree D = sz3-1
% and P0 = A(:,:,end),..., PD = A(:,:,1). 
%
% Those who prefer the coefficients to be ordered according to ascending 
% powers should use POL.

%   The internal structure of the polynomial matrix object P is as follows:
%
%   P.d is the degree of P(s);
%
%   P.s is the size of P(s), a two-element row vector [nrow, ncol] containing
%       the number of rows and columns;
%
%   P.c is a three-dimensional array containing matrix coefficients of P(s).
%       More specifically, P.c(:,:,1) = P0, P.c(:,:,2) = P1 and
%       P.c(:, :, P.d+1) = Pd.
%
%   P.v is the variable string;
%
%   P.u is a user data field;

%       Author(s):  S. Pejchova 05-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 11-Oct-1999 12:00:00  J.Jezek  $
%       $Revision: 3.0 $  $Date: 08-Nov-1999 14:45:11  S.Pejchova  $
%                         $Date: 31-May-2000 16:00:00  J.Jezek  $
%                         $Date: 03-Oct-2002  S.Pejchova  $

na = nargin;
if na & isa(varargin{1},'pol'),

   % Quick exit for POL(P) with P of class POL
   P = varargin{1};
   if na > 1,
      error('Use PPROPS to modify the properties of the POL object.');
   end;
   return
   
elseif ~na,
   P = pol; return;
   
else
   Po = varargin{1}; dx=0;
   try, 
       if ndims(Po)==3, P = pol(Po); 
       else, P = pol(Po,varargin{2:na}); 
       end; 
   catch,
       fx=findstr(lower(lasterr),'degree');
       if ~isempty(fx), 
           dx =[];
           try, P = pol(Po); dx=0; catch,  end; 
       end;
       if isempty(dx), error(peel(lasterr)); end;
   end;
   dd = P.d; Pc = P.c;
   if na>1, 
       dx = varargin{2};
       if ~isa(dx,'double') | length(dx)~=1 | ~isfinite(dx), 
           dx = 0; 
       end; 
       if dx>0,
            Pc(:,:,dd+2:dx+1) = 0; 
       end; 
   end;
   if ndims(Pc)==3, 
       Pc(:,:,:) = Pc(:,:,end:-1:1); 
       try, P = pol(Pc,varargin{2:na}); 
       catch, error(peel(lasterr));
       end; 
   end; 
end;
   
%end .. lop
