function C = rchar(P,n)
%RCHAR  Convert a polynomial object to a string in a reverse order of degrees
%
% It is dual to the function POL/CHAR.
% The commmand
%    RCHAR(P) 
% is the default. It looks like RCHAR(P,4) but uses roughly four digits of 
% precision and an exponent if required. The commmand
%    RCHAR(P,N) 
% converts the scalar polynomial P into a string. If P is a polynomial 
% matrix then the result is a cell array of the same size consisting of 
% strings that correspond to the polynomial elements of P. Each scalar 
% coefficient is taken with the precision of at most N digits. The 
% commmand
%    RCHAR(P,FORMAT) 
% works like RCHAR(P,N) but uses the format string FORMAT for each scalar
% coefficient (see SPRINTF for details). The commmand
%    RCHAR(P,'RAT') 
% is a special case of RCHAR(P,FORMAT) and uses rational aproximation of 
% the coefficients.

%      Author(s): S. Pejchova 24-6-99
%      Copyright (c) 1999 by Polyx, Ltd.
%      $Revision: 3.0 $  $Date: 24-June-1999 10:51:43   $
%                        $ 19-Aug-2001   J.Jezek  arg checking  $
% Effect on other properties:
% C is a standard Matlab cell array of strings.

ni = nargin;
if ni < 1,
   error('Not enough input arguments.');
end;
eval('P = pol(P);', 'error(peel(lasterr));');
[s1,s2]=size(P); Pc=P.c; v=P.v;
C=cell(s1,s2);
if ~isempty(C),
   if isinf(P.d), Pc=zeros(s1,s2);
   end;
   d0=size(Pc,3)-1;
   for i=1:s1, for j=1:s2,
     ePc = Pc(i,j,d0+1:-1:1);  ePc = ePc(:,:);
     if any(ePc),
        s=''; d=d0;
        for x=ePc,
          x1=0; sxx=''; sx1 = ''; sx2 = '';
          if x~=0,
            if ~isempty(s),    
              if ((real(x)>=0)|(isnan(x))), s = [s ' + '];
              else, s = [s ' - ']; x = -x; x1=1;
              end;
            elseif ((real(x)>=0)|(isnan(x))), s=' ';
            else, s='-'; x=-x; x1=1;
            end;
            if x~=1|d==0,
              if ni==1, sxx = num2str(x);
              elseif ischar(n) & strcmp(n,'rat'),
                 [n1,d1]=rat(x);
                 sxx = num2str(n1);
                 if abs(d1)~=1, sx1 = ['/',int2str(d1)];
                 end;
                 if d>0,sx2 = '*';
                 end;   
              elseif (ischar(n) & ~isempty(n)) | ...
                     (isnumeric(n) & length(n)==1),
                 eval('sxx = num2str(x,n);','error(peelf(lasterr));');
              else
                 error('Invalid 2nd argument.');
              end;
              if ~isreal(x)&(d>0|~isempty(sx1)|x1),
                 sxx = ['(',sxx,')'];
              end;
              if ~strcmp([sxx sx1],'1')|(d==0),
                 s = [s sxx sx1 sx2];
              end;
            end; %if x~=1|d==0
            if d>=2,
              if strcmp(v,'z^-1'), s=[s, 'z^-', int2str(d)];
              else, s=[s, v, '^', int2str(d)];
              end;
            elseif d==1, s=[s v];
            end;
          end; %if x~=0
          d=d-1;
        end; % for x=ePc
     else,
        s=' 0';
     end; % any(ePc)
     C{i,j}=s;
   end; end;
end;
if (s1==1 & s2==1), C=C{1,1}; end;

%end .. rchar
