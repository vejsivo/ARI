function C = char(T,n)
%CHAR  Convert a two-sided polynomial object to a string
%
% The commmand
%    CHAR(T) 
% converts the scalar polynomial T into a string. If T is a two-sided polynomial 
% matrix then the result is a cell array of the same size consisting of 
% strings that correspond to the polynomial elements of T. Each scalar 
% coefficient is taken with respect to the basic Matlab format.
% The commmand
%    CHAR(T,N) 
% converts each scalar coefficient into a string representation with 
% a maximum N digits of precision. The commmand
%    CHAR(T,FORMAT) 
% works like CHAR(T,N) but uses the format string FORMAT for each scalar
% coefficient (see SPRINTF for details). The commmand
%    CHAR(T,'RAT') 
% is a special case of CHAR(T,FORMAT) and uses rational aproximation of 
% the coefficients.
%    CHAR(T,'ROOTC') 
% uses zero/pole/gain representation of polynomial entries 
%    CHAR(T,'ROOTR') 
% is a pretty-form of CHAR(T,'ROOTC') without complex numbers (the complex 
% pairs are multiplied).

%      Author(s): S. Pejchova  08-7-99
%      Copyright (c) 1998 by Polyx, Ltd.
%      $Revision: 3.0 $  $Date: 17-Sep-1999 11:56:43   $
%      $Revision: 3.0 $  $Date: 05-Jun-2000 11:18:43   $

% Effect on other properties:
% C is a standard Matlab cell array of strings.

global PGLOBAL;

ni = nargin;
if ni == 1,
   n = [];
elseif ~(isa(n,'double') | isa(n,'char')),
   error('Invalid 2nd argument.');
end;

s1 = T.s(1); s2 = T.s(2);
Tp = T.p; Tc = Tp.c; offst = T.o; Tr = T.r;

orr = 0;
if strcmp(PGLOBAL.ORDER,'reverse'),
   orr = 1;
end;

C = cell(s1,s2);
if ~isempty(C),
   if isinf(Tr),
      Tc = zeros(s1,s2);
   end;
   for i = 1:s1, for j = 1:s2,
     eTc = Tc(i,j,:);  eTc = eTc(:,:);
     if any(eTc),
        s = ''; d0 = offst; eTc = eTc(1:max(find(eTc)));
        if strncmp(n,'root',4) & length(find(eTc))>1,
           r1 = roots(fliplr(eTc));
           r1 = sort(r1); [f_x,ix] = sort(real(r1));
           r1 = r1(ix); s = ' '; eTc_e = eTc(end);
           if isreal(eTc_e), 
              if (eTc_e < 0),
                 s = '-';
              end; 
              if abs(eTc_e)~=1,
                 s = [s,fn2str(abs(eTc_e))];
              end;
           else,
              s = [s,'(',fn2str(eTc_e),')']; 
           end;
           fr_z = find(r1==0);
           r1(fr_z) = []; l_fr_z = length(fr_z)+d0;
           if l_fr_z,
              s = [s,'z'];
              if l_fr_z~=1,
                 s = [s,'^',int2str(l_fr_z)];
              end;
           end;
           lr1 = length(r1);
           if lr1>0,
              ixc = 1; r1 = -r1;
              while ixc<=lr1,
                 if strcmp(n,'rootr') & ixc<lr1 & (~isreal(r1(ixc))) &...
                       (~isreal(r1(ixc+1))) & (abs(r1(ixc)-conj(r1(ixc+1)))<eps),
                    sgr1 = [fn2str(abs(2*real(r1(ixc)))),'z'];
                    if real(r1(ixc))>eps,
                       sgr1 = ['+',sgr1];  
                    elseif real(r1(ixc))<-eps,
                       sgr1 = ['-',sgr1];
                    else,
                       sgr1 = ''; 
                    end;
                    s = [s,'(z^2',sgr1,'+',...
                          fn2str((real(r1(ixc)))^2+(imag(r1(ixc)))^2),')'];
                    ixc = ixc+2;
                 else,
                    if abs(real(r1(ixc)))<eps,
                       sgr1 = [fn2str(imag(r1(ixc))),'i'];
                    else,
                       sgr1 = fn2str(r1(ixc));
                    end;
                    if ~strcmp(sgr1(1),'-'),
                       sgr1 = ['+',sgr1];
                    end;
                    s = [s,'(z',sgr1,')'];
                    ixc = ixc+1;
                 end;
              end;
           end;
        else,
           if orr,
              d0 = offst+length(eTc)-1; eTc = fliplr(eTc);
              eTc = eTc(1:max(find(eTc)));
           end;
           for x = eTc,
              d = abs(d0); sxx = '';
              if x~=0,
                 if ~isempty(s),
                    if isreal(x) & (x<0) & (~isnan(x)),
                       s = [s ' - '];  x = -x;
                    else,
                       s = [s ' + '];
                    end;
                 elseif x==-1,
                    s = '-'; x = -x;
                 elseif isnan(real(x)) | (real(x)>=0) | (~isreal(x)&d>0),
                    s = ' ';
                 end;
                 if x~=1 | d==0,
                    if ni==1 | strncmp(n,'root',4), 
                       sxx = fn2str(x);
                    elseif strcmp(n,'rat'),
                       [n1,d1] = rat(x); sxx = num2str(n1);
                       if abs(d1)~=1, 
                          if ~isreal(n1),
                             sxx = ['(',sxx,')'];
                          end;
                          sxx = [sxx,'/',int2str(d1)];
                       end;
                    else,
                       sxx = num2str(x,n);
                    end;
                    if ~isreal(x) & (d>0),
                       sxx = ['(',sxx,')'];
                    end;
                    if ~strcmp(sxx,'1') | (d==0),
                       s = [s sxx];
                    end;
                 end; %if x~=1|d==0
                 if d,
                    s = [s,'z'];
                 end;
                 if d>=2 | d0<0 ,
                    s = [s, '^'];
                    if d0<0,
                       s = [s, '-'];
                    end;
                    s = [s, int2str(d)];
                 end;
                 if length(s)>1,
                    if strcmp(s(1:2),'(-'),
                       s = [' ',s];
                    end; 
                 end;
              end; %if x~=0
              if orr,
                 d0 = d0-1;
              else
                 d0 = d0+1;
              end;
           end; % for x=eTc
        end; %if strncmp(n,'root',4)
     else,
        s = ' 0';
     end; % any(eTc)
     C{i,j} = s;
   end; end;
end;
if (s1==1 & s2==1),
   C = C{1,1};
end;

%Subfunction - Convert number to string with respect to the basic Matlab format
function s_f = fn2str(nn)
fr_b=get(0,'Format');
s_f='';
switch fr_b,
case 'short',
   if (nn==round(nn)& abs(nn)<1e10),
      s_f=num2str(nn,'%11.0f');
   elseif (nn~=round(nn)& (abs(nn)<1000) & (abs(nn)>0.001)),
      s_f=num2str(nn,'%11.4f');
   else,
      s_f=num2str(nn,'%11.4e');
   end;
case {'shortE','hex'},
   if (nn==round(nn)& abs(nn)<1e10),
      s_f=num2str(nn,'%11.0f');
   else,
      s_f=num2str(nn,'%11.4e');
   end;
case 'shortG',
   s_f=num2str(nn,'%11.5g');
case 'long',
   if (nn==round(nn)& abs(nn)<1e10),
      s_f=num2str(nn,'%21.0f');
   elseif (nn~=round(nn)& (abs(nn)<100) & (abs(nn)>0.001)),
      s_f=num2str(nn,'%21.15f');
   else,
      s_f=num2str(nn,'%21.15e');
   end;  
case 'longE',
   if (nn==round(nn)& abs(nn)<1e15),
      s_f=num2str(nn,'%21.0f');
   else,
      s_f=num2str(nn,'%21.15e');
   end;
case 'longG',
   s_f=num2str(nn,'%21.15g');
case 'bank',
   s_f=num2str(nn,'%11.2f');
case 'rational',
   [n1,d1]=rat(nn);
   sd1='';
   sn1 = num2str(n1);
   if abs(d1)~=1, 
      sd1 = ['/',int2str(d1)];
      if ~isreal(n1), sn1=['(',sn1,')']; end;
   end;
   s_f=[sn1,sd1];
case '+',
   if nn>0, s_f='+';
   elseif nn<0, s_f='-';
   end;
end;
   

%end .. @tsp/char
