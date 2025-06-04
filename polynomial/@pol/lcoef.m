function [L,D] = lcoef(A,varargin)
%LCOEF  Leading coefficient of polynomial
%
% L = LCOEF(A)       default, the same as LCOEF(A,'mat')
% L = LCOEF(A,'mat') returns the leading coefficient matrix of A
% L = LCOEF(A,'ent') returns the scalar leading coefficients of
%                    the entries of A
% L = LCOEF(A,'row') returns the row leading coefficient matrix of A
% L = LCOEF(A,'col') returns the column leading coefficient matrix of A
% L = LCOEF(A,'dia') for para-Hermitian polynomial matrix A,
%                     returns the diagonally leading coefficient matrix
%
% For polynomials in 'z' or 'z^-1', it is possible to specify
% (by an optional input argument VAR) that the leading coefficient is to be
% understood by the highest power of 'z' or 'z^-1'.
%
% The second output argument D in all cases returns the corresponding
% matrix or vector of degrees. D is the same as the first output argument
% of the function DEG.
%
% See also POL/DEG.

%       Author(s):  S. Pejchova, M. Sebek 17-3-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 22-Apr-1999 12:20:34   $
%       $Revision: 3.0 $  $Date: 13-Oct-1999 12:00:00   J. Jezek  $
%                         $Date: 13-Oct-2002            J. Jezek  $
%                         $Date: 17-Nov-2002  z,z^-1    J. Jezek  $

% Effect on other properties:
% L and D are standard Matlab matrices.

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
      eval('L = tcoef(A,string);','error(peel(lasterr));');
   else
      eval('[L,D] = tcoef(A,string);','error(peel(lasterr));');
   end;
   D = -D;
   return;
end      

Ac = A.c; [As1,As2,Ad] = size(A);
L = zeros(As1,As2);

if ~isempty(Ad) & Ad>0,
   a = repmat(reshape(1:Ad+1,[1 1 Ad+1]),[As1 As2]);
end;

switch string
 case 'mat',
      if isempty(Ad), D =[]; return; end;
      if ~isinf(Ad), L = Ac(:,:,Ad+1); end;
      D = Ad;
   
 case 'ent',
      if isempty(Ad), D = zeros(As1,As2); return; end; 
      if ~isinf(Ad),       
%        Do = Ac&Ac; L = Ac;   does not work for Ac complex
%                              correction J.Jezek 13-Oct-2002
         Do = Ac~=0; L = Ac;
         if Ad>0,  Do = Do.*a; Do = max(Do,[],3); end;
         D = Do-1;
         D(D<0) = -Inf;
         if Ad > 0,
            Do = ((repmat(Do,[1,1,Ad+1]))==a);
            L = sum((Do.*Ac),3);
         end;
      elseif no==2,
         D = repmat(-Inf,[As1,As2]);
      end;

case 'row',
      if isempty(Ad), D = zeros(As1,0); return; end;
      if ~isinf(Ad),
%        Do = Ac&Ac; L = Ac;
         Do = Ac~=0; L = Ac;
         if Ad>0, Do = Do.*a;  end;
         Do = max(Do(:,:),[],2);  D = Do-1;
         D(D<0) = -Inf;
         if Ad > 0,
            Do = ((repmat(Do,[1,As2,Ad+1]))==a);
            L = sum((Do.*Ac),3);
         end;        
      elseif no==2,
         D = repmat(-Inf,[As1,1]);
      end;

 case 'col',
      if isempty(Ad), D = zeros(0,As2); return; end;
      if ~isinf(Ad),
%        Do = Ac&Ac; L = Ac;
         Do = Ac~=0; L = Ac;
         if Ad>0, Do = Do.*a; Do = max(Do,[],3); end;
         Do = max(Do,[],1);  D = Do-1;
         D(D<0)=-Inf;            
         if Ad > 0,
            Do = ((repmat(Do,[As1,1,Ad+1]))==a);
            L = sum((Do.*Ac),3);
         end;
      elseif no==2,
         D = repmat(-Inf,[1,As2]);
      end;
 case 'dia',
      global PGLOBAL;
      if As1~=As2,
         error('Matrix is not square.');
      end;
      if isempty(Ad), D =[]; return; end;
      if norm(A-A','blk',inf) > PGLOBAL.ZEROING,
         error('Matrix is not para-Hermitian.');
      end;

      if ~isinf(Ad),
%        Do=Ac&Ac;
         Do = Ac~=0;
         if Ad>0, Do = Do.*a; Do = max(Do,[],3); end;
         D = (Do-1).*((Do-1)>=0);
         halfdiag = 1;
         for i = 1:As1, for j = i+1:As1,
             if 2*D(i,j) > D(i,i)+D(j,j), halfdiag = 0; end;
         end; end;
	
         if halfdiag,
            % if halfdiag = 1 diagonal degrees are half the
            % degrees of diagonal entries

            D = diag(D)/2;

         else
            warning('The half diagonal degrees are not unique.');
            % computation of the smallest diagonal degrees sdiag(i)
            % such that D(i,j) <= di+dj
            ddiag = [D(1,1)/2 zeros(1,As1-1)];
            for i = 2:As1,
                ddiag(i) = max([D(i,1:i-1)-ddiag(1:i-1) D(i,i)/2]);
            end;
            D = ddiag';
         end; % if halfdiag

         for i = 1:As1, for j = 1:As1,
             if D(i)+D(j) <= Ad,
                b = Ac(i,j,(D(i)+D(j)+1));
                L(i,j) = (-1)^D(i)*b;
             end;
         end; end;

      elseif no==2,
         D = repmat(-Inf,[As1,1]);
      end;

 otherwise,
   error('Invalid command option.');

end; %switch string

%end .. @pol/lcoef
