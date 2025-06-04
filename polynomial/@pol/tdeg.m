function [D,T] = tdeg(A,varargin)
%TDEG    Trailing degree of polynomial
%
% D = TDEG(A)       default, the same as TDEG(A,'mat').
% D = TDEG(A,'mat') returns the trailing degree of polynomial matrix A.
% D = TDEG(A,'ent') returns the matrix of trailing degrees of the A entries.
% D = TDEG(A,'row') returns the column vector of row trailing degrees of A.
% D = TDEG(A,'col') returns the row vector of column trailing degrees of A.
% D = TDEG(A,'dia') for para-Hermitian polynomial matrix A,
%                    returns the vector of half diagonal trailing degrees of A.
%
% For polynomials in 'z' or 'z^-1', it is possible to specify
% (by an optional input argument VAR) that the trailing degree is to be
% understood by the lowest power of 'z' or 'z^-1'.
%
% The second output argument T in all cases returns corresponding
% matrix or vector of trailing coefficients, the same as first output argument
% in function TCOEF.
%
% See also POL/TCOEF.

%          Author: J. Jezek  13-10-99
%          Copyright (c) 1999 by Polyx, Ltd.
%          $ Revision 3.0 $   $ Date 06-Apr-2000 $
%                             $ Date 20-Jul-2000 $
%                             $ Date 07-Nov-2000 $
%                             $ Date 17-Nov-2002  z,z^-1  $
% Effect on other properties:
% D and T are standard Matlab matrices.

if ~isa(A,'pol'),
   error('Some argument but not 1st is invalidly pol.');
end;
Av = A.v; recip = logical(0);

string = 'mat';
li = length(varargin);
if li>0,
   for i = 1:li,
      arg = varargin{i};
      if isa(arg,'pol'),
         [vs1,vs2,vd] = size(arg);
         if all([vs1,vs2,vd]==1) & all(arg.c(:,:)==[0 1]),
            arg = arg.v;
         else
            error(['Invalid ',nth(i+1),' argument.']);
         end;
      end;
      if ischar(arg),
         if strcmp(arg,'zi'), arg = 'z^-1';
         end;
         I = strmatch(arg,{'s';'p';'z^-1';'d';'z';'q'},'exact');
         if ~isempty(I),
            if ~isempty(Av) & ~strcmp(Av,arg),
               if (strcmp(Av,'z') & strcmp(arg,'z^-1')) | ...
                     (strcmp(Av,'z^-1') & strcmp(arg,'z')),
                  recip = logical(1);
               else
                  error('Invalid variable symbol.');
               end;
            end; 
         else
            string = arg;
         end;
      else
         error(['Invalid ',nth(i+1),' argument.']);
      end;
   end;
end;        
   
no = nargout;
if recip,
   if no<2,
      eval('D = deg(A,string);','error(peel(lasterr));');
   else
      eval('[D,T] = deg(A,string);','error(peel(lasterr));');
   end;
   D = -D;
   return;
end;      

if ~isempty(A),
   Ad = A.d;
   if strcmp(string,'dia') & ~(floor(Ad/2)*2==Ad), Add = Ad+1;
   else Add = Ad;
   end;
   A = rev(A,Add);
end;

D = 0; T = 0;
eval('[D,T] = deg(A,string);','error(peel(lasterr));');
Ds1 = size(D,1); Ds2 = size(D,2);
for i = 1:Ds1,
   for j = 1:Ds2,
      if isfinite(D(i,j)),
         if strcmp(string,'dia'),
            D(i,j) = floor(Add-D(i,j))/2;
         else
            D(i,j) = Add - D(i,j);
         end;
      else
         D(i,j) = inf;
      end;
   end;
end;

%end .. @pol/tdeg

   