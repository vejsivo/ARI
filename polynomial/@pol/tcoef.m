function [T,D] = tcoef(A,varargin)
%TCOEF   Trailing coeficient of polynomial
%
% T = TCOEF(A)         default, the same as TCOEF(A,'mat')
% T = TCOEF(A,'mat')   returns the trailing coefficient matrix of A
% T = TCOEF(A,'ent')   returns the scalar trailing coefficients
%                        of the entries of A
% T = TCOEF(A,'row')   returns the row trailing coefficient matrix of A
% T = TCOEF(A,'col')   returns the column trailing coefficient matrix of A
% T = TCOEF(A,'dia')   returns for a para-Hermitian polynomial matrix A
%                        the diagonal trailing coefficient matrix
%
% For polynomials in 'z' or 'z^-1', it is possible to specify
% (by an optional input argument VAR) that the trailing coefficient is to be
% understood by the highest power of 'z' or 'z^-1'.
%
% The second output argument D in all cases returns the corresponding
% matrix or vector of trailing degrees. D is the same as the first output
% argument of the function TDEG.
%
% See also POL/TDEG.

%          Author: J. Jezek  13-10-99
%          Copyright (c) 1999 by Polyx, Ltd.
%          $ Revisioin $  $ Date 07-Nov-2000 $
%                         $ Date 17-Nov-2002 $

% Effect on other properties:
% T and D are standard Matlab matrices.

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
      eval('T = lcoef(A,string);','error(peel(lasterr));');
   else
      eval('[T,D] = lcoef(A,string);','error(peel(lasterr));');
   end;
   D = -D;
   return;
end;      

Ad = A.d;
Ac = A.c; Ac = flipdim(Ac,3); A.c = Ac;
if ~isempty(Ad) & isfinite(Ad),
   A = pol(Ac(:,:),Ad);
end;
T = 0; D = 0;
eval('[T,D] = lcoef(A,string);','error(peel(lasterr));');
Ds1 = size(D,1); Ds2 = size(D,2);
for i = 1:Ds1,
   for j = 1:Ds2,
      if isfinite(D(i,j)),
         D(i,j) = Ad-D(i,j);
      else
         D(i,j) = inf;
      end;
   end;
end;

%end .. @pol/tcoef

   