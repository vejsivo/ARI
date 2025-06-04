function C = char(P,n)
%CHAR  Convert a polynomial object to a string
%
% The commmand
%    CHAR(P) 
% converts the scalar polynomial P into a string. If P is a polynomial 
% matrix then the result is a cell array of the same size consisting of 
% strings that correspond to the polynomial elements of P. Each scalar 
% coefficient is taken with respect to the basic Matlab format.
% The commmand
%    CHAR(P,N) 
% converts each scalar coefficient into a string representation with 
% a maximum N digits of precision. The commmand
%    CHAR(P,FORMAT) 
% works like CHAR(P,N) but uses the format string FORMAT for each scalar
% coefficient (see SPRINTF for details). The commmand
%    CHAR(P,'RAT') 
% is a special case of CHAR(P,FORMAT) and uses rational aproximation of 
% the coefficients.
%    CHAR(P,'ROOTC') 
% uses zero/pole/gain representation of polynomial entries 
%    CHAR(P,'ROOTR') 
% is a pretty-form of CHAR(P,'ROOTC') without complex numbers (the complex 
% pairs are multiplied).

%      Author(s): S. Pejchova, M. Sebek 03-4-98
%      Copyright (c) 1998 by Polyx, Ltd.
%      $Revision: 3.0 $  $Date: 05-Jun-2000 11:09:43   $
%                     $  $Date: 13-Mar-2003 S.Pejchova    $

% Effect on other properties:
% C is a standard Matlab cell array of strings.

global PGLOBAL;

ni = nargin;
if ni == 1,
   n = [];
elseif ~(isa(n,'double') | isa(n,'char')),
   error('Invalid 2nd argument.');
end;

s1 = P.s(1); s2 = P.s(2);
Pc = P.c; v = P.v; v1 = v;
if strcmp(v1,'z^-1'), 
   v1 = 'z^-';
else,
   v1 = [v1,'^'];
end;

orr = 0;
if strcmp(PGLOBAL.ORDER,'reverse'),
   orr = 1;
end;

C = cell(s1,s2);
if ~isempty(C),
   if isinf(P.d),
      Pc = zeros(s1,s2);
   end;
   for ii = 1:s1, for jj = 1:s2,
     ePc = Pc(ii,jj,:);  ePc = ePc(:).'; 
     if any(ePc),
        s = ''; d = 0; ePc = ePc(1:max(find(ePc)));
        if strncmp(n,'root',4) & length(ePc)>1,
           r1 = roots(fliplr(ePc));
           r1 = sort(r1); [f_x,ix] = sort(real(r1));
           r1 = r1(ix); s = ' '; ePc_e = ePc(end);
           if isreal(ePc_e), 
              if (ePc_e < 0),
                 s = '-';
              end;
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % if abs(ePc_e)~=1,
              if abs(abs(ePc_e)-1)>eps,
                 s = [s,fn2str(abs(ePc_e))];
              end;
           else,
              s = [s,'(',fn2str(ePc_e),')']; 
           end;
           fr_z = find(r1==0);
           if length(fr_z)>1, 
              s = [s,v1,int2str(length(fr_z))];
              r1(fr_z) = [];
           elseif length(fr_z)==1,
              s = [s,v];
              r1(fr_z) = [];
           end;
           lr1 = length(r1);
           if lr1>0,
              ixc = 1; r1 = -r1;
              while ixc<=lr1,
                 if strcmp(n,'rootr') & ixc<lr1 & (~isreal(r1(ixc))) &...
                       (~isreal(r1(ixc+1))) & (abs(r1(ixc)-conj(r1(ixc+1)))<eps),
                    sgr1 = [fn2str(abs(2*real(r1(ixc)))),v];
                    if real(r1(ixc))>eps,
                       sgr1 = ['+',sgr1];  
                    elseif real(r1(ixc))<-eps,
                       sgr1 = ['-',sgr1];
                    else,
                       sgr1 = ''; 
                    end;
                    s = [s,'(',v1,'2',sgr1,'+',...
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
                    s = [s,'(',v,sgr1,')'];
                    ixc = ixc+1;
                 end;
              end;
           end;
        else,
           if orr,
              d = length(ePc)-1; ePc = fliplr(ePc);
              ePc = ePc(1:max(find(ePc)));
           end;
           for x = ePc,
              sxx = '';
              if x~=0,
                 if ~isempty(s),
                    if isreal(x) & (x<0) & (~isnan(x)),
                       s = [s ' - '];  x = -x;
                    else,
                       s = [s ' + '];
                    end;
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 % elseif x==-1,
                 elseif abs(x+1)<= eps,
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
                 if d>=2,
                    s = [s, v1, int2str(d)];
                 elseif d==1,
                    s = [s v];
                 end;
                 if length(s)>1,
                    if strcmp(s(1:2),'(-'),
                       s = [' ',s];
                    end; 
                 end;
              end; %if x~=0
              if orr,
                 d = d-1;
              else
                 d = d+1;
              end;
           end; % for x=ePc
        end; %if strncmp(n,'root',4)
     else,
        s = ' 0';
     end; % any(ePc)
     C{ii,jj} = s;
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
end

%end .. @pol/char
