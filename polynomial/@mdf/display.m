function display(R)
%DISPLAY  Display matrix-denominator fraction
%
% Syntax
%    DISPLAY(R)
% The display format is set by macro PFORMAT.

%       Author(s):  S. Pejchova  25-01-00
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 27-Jan-2000 15:36:27   $
%                         $Date: 20-Mar-2000 12:00:00 J.Jezek $
%                         $Date: 11-Jul-2000 10:38:00 S.Pejchova $
%                         $Date: 14-Sep-2000  S.Pejchova  $
%                         $Date: 14-Oct-2002  J.Jezek     $
%                         $Date: 31-Mar-2003  J.Jezek     $

global PGLOBAL

row=0;  fr=[];  ls=1;  str_crp='';

eval('fr=PGLOBAL.FORMAT;',...
   'error(''Use PINIT to initialize the Polynomial Toolbox.'');');
if strcmp(get(0,'FormatSpacing'),'compact'), ls=0;
end;
fr_b=get(0,'Format');
vr=R.frac.v;  s1=R.frac.s(1);  s2=R.frac.s(2);
iname=inputname(1);
if strcmp(fr,'nice'),
   [aR1,cR1]=declass(R);
   if ~strcmp(cR1,'mdf'),
      eval([iname,'=aR1;  display(',iname,');']);
      return; 
   end; 
end;

iname=[iname,' = '];
if ls, iname=char('   ',iname);
end;

if ~strcmp(R.frac.c,'cop?'),
   str_crp='coprime';
   if strcmp(R.frac.c,'ncop'),
      str_crp=['non',str_crp];
   end;
end;
if ~strcmp(R.frac.r,'red?'),
   s_r='reduced';
   if strcmp(R.frac.r,'nred'),
      s_r=['non',s_r];
   end;
   if isempty(str_crp), str_crp=s_r; 
   else, str_crp=[str_crp,', ',s_r]; 
   end;
end;
if ~strcmp(R.frac.p,'prop?'),
   s_r='proper';
   if strcmp(R.frac.p,'nprop'),
      s_r=['im',s_r];
   end;
   if isempty(str_crp), str_crp=s_r; 
   else, str_crp=[str_crp,', ',s_r]; 
   end;
end;
if strcmp(fr,'nice'), str_crp='';
end;
flags1=1;
if isempty(str_crp), flags1=0;
end;

if isempty(R.frac.num.d),
   S_out=['   Empty matrix-denominator fraction (numerators==denominators): ',...
         int2str(s1),'-by-',int2str(s2)];
   if flags1 ,S_out=[S_out,'   ',str_crp];
   end;
   if ls, 
      disp(char(iname,' ',S_out,' ')); 
   else,
      disp(char(iname,S_out));
   end;
   return;
end;

in_coef=strmatch(fr,['coef ';'rcoef';'block'],'exact');
if ~isempty(in_coef), fr='symb';
end;


n_lg=[];
switch fr,
case {'symbs','nice'}, n_lg=2;  
case 'symbr', n_lg='rat';
case {'rootr','rootc'}, n_lg=fr;
end;
if isempty(n_lg), 
   Cn=char(R.frac.num);
   Cd=char(R.frac.den); 
else, 
   Cn=char(R.frac.num,n_lg);
   Cd=char(R.frac.den,n_lg);
end;
C = cell(s1,s2);
if (s1==1&s2==1), Cn=cellstr(Cn);  Cd=cellstr(Cd);
end; 
for ii=1:s1,
   for jj=1:s2,
      Cdij=Cd{ii,jj}; Cnij=Cn{ii,jj};
      if strcmp(Cdij(1),' '), Cdij=Cdij(2:end);
      end;
      if strcmp(Cnij(1),' '), Cnij=Cnij(2:end);
      end;
      lCd=length(Cdij); lCn=length(Cnij);
      if ~strcmp(Cdij,'1'),
         C{ii,jj}=strvcat([repmat(' ',[1,max(0,ceil((lCd-lCn)/2))]),Cnij],...
            repmat('-',[1,max(lCn,lCd)]),...
            [repmat(' ',[1,max(0,ceil((lCn-lCd)/2))]),Cdij],' ');
      else
         C{ii,jj}=strvcat(' ',Cnij,' ',' ');
      end;  
   end; 
end;

const1=0;
c1=repmat([' '],[4*s1,4]);
if (~R.frac.num.d)&(~R.frac.den.d),
   iname1=['Constant matrix-denominator fraction: ',int2str(s1),'-by-',int2str(s2)];
   if flags1, iname1=[iname1,'   ',str_crp]; flags1=0; end;
   iname1=char(iname1,iname); iname=iname1;
   if ls, iname=char(' ',iname); end;
   const1=1;
else,
   iname1=['Matrix-denominator fraction in ',vr,': ',int2str(s1),...
         '-by-',int2str(s2)];
   if flags1, iname1=[iname1,'   ',str_crp]; end;
   iname1=char(iname1,iname);
end;
S_out=[];
if ls,  iname1=char(' ',iname1); end;
if strcmp(fr,'nice'), iname1=iname; end;
if s2==1,
   S_out=[c1,char(C{:,1})]; 
   if ls, S_out=char(' ',S_out,' '); end;
   S_out=char(iname,S_out);
else,
   S=char(C{:,1}); S=[c1,S];
   col=1; sumlg=size(S,2);
   for i=2:s2+1,
      if i<=s2,
         Saux=char(C{:,i}); Saux=[c1,c1(:,1:2),Saux];
         onelg=size(Saux,2);
      end;
      if (sumlg+onelg < 79)&(i<=s2),
         col=[col,i]; S=[S,Saux]; sumlg=sumlg+onelg;
      else,
         if length(col)~=s2,
            if any(col==1),
               S_out=iname1;  flags1=0; 
               if ls, S_out=char(S_out,' '); end;
            end;
            lgcol=length(col);
            if lgcol==1,
               S_col=['   Column ',int2str(col)];
            else,
               S_col=['   Columns ',int2str(col(1)),' through ',int2str(col(lgcol))];
            end;
            if ls, S_col=char(S_col,' '); end;
            S_out=char(S_out,S_col);
         else,
            S_out=iname;
            if ls, S_out=char(S_out,' '); end;
         end;
         if isempty(S_out)&~ls,
            S_out=S;
         else,
            S_out=char(S_out,S); 
         end;
         if ls, S_out=char(S_out,' '); end;
         col=i; S=Saux(:,3:end); sumlg=onelg;
      end;
   end; % for i=1:s2
end; % if (s2==1...
disp(S_out(1:end-1,:));
if flags1, disp(['    ',str_crp]); end;
if ls,  disp(' '); end;

%end .. @mdf/display
