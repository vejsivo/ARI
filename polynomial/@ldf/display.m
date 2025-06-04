function display(F)
%DISPLAY   Display left-denominator fraction
%
% Syntax
%    DISPLAY(F)
% The display format is set by macro PFORMAT.

%       Author(s):  S. Pejchova  05-01-00
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 16-Jan-2000 12:00:27   $
%                         $Date: 21-Apr-2000 12:00:00 J.Jezek $
%                         $Date: 14-Sep-2000  S.Pejchova  $
%                         $Date: 14-Oct-2002  J.Jezek     $


global PGLOBAL

row=0;  ls=1;  str_crp='';

eval('fr=PGLOBAL.FORMAT;',...
   'error(''Use PINIT to initialize the Polynomial Toolbox.'');');
if strcmp(get(0,'FormatSpacing'),'compact'), ls=0; end;
iname=inputname(1);
if strcmp(fr,'nice'),
   [aF1,cF1]=declass(F);
   if ~strcmp(cF1,'ldf'),
      eval([iname,'=aF1;  display(',iname,');']);
      return; 
   end; 
end;

if ~strcmp(F.frac.c,'cop?'),
   str_crp='coprime';
   if strcmp(F.frac.c,'ncop'),
      str_crp=['non',str_crp];
   end;
end;
if ~strcmp(F.frac.r,'red?'),
   s_r='reduced';
   if strcmp(F.frac.r,'nred'),
      s_r=['non',s_r];
   end;
   if isempty(str_crp), str_crp=s_r; 
   else, str_crp=[str_crp,', ',s_r]; 
   end;
end;
if ~strcmp(F.frac.p,'prop?'),
   s_r='proper';
   if strcmp(F.frac.p,'nprop'),
      s_r=['im',s_r];
   end;
   if isempty(str_crp), str_crp=s_r; 
   else, str_crp=[str_crp,', ',s_r]; 
   end;
end;
if strcmp(fr,'nice'), str_crp=''; end;

str_L=pdisp(F.frac.den);
str_R=pdisp(F.frac.num);
fe1=findstr('Empt',str_L(ls+1,:));
if isempty(fe1), 
   fe1=findstr('matr',str_L(ls+1,:));
   if ~isempty(fe1), str_L(1:ls+1,:)=[]; end; 
end;
fe1=findstr('Empt',str_R(ls+1,:));
if isempty(fe1), 
   fe1=findstr('matr',str_R(ls+1,:));
   if ~isempty(fe1), str_R(1:ls+1,:)=[]; end; 
end;
str_L=str_L(:,1:max(find(any(str_L~=' ',1))));
str_R=str_R(:,1:max(find(any(str_R~=' ',1))));

[l1,l2]=size(str_L);
[r1,r2]=size(str_R);
if (l1<21) & (r1<21) & (l2+r2<72), row=1; end;

if row,
   s1=max(l1,r1); zn=s1;
   if ls, zn=max(zn-2,1); end;
   zn=min(5,zn);
   s2=l2+r2+6+zn;
   S_out=repmat(' ',[s1,s2]);
   w1=fix((s1-l1)/2); w2=fix((s1-r1)/2); w3=fix((s1-zn)/2);
   S_out(w1+1:w1+l1,1:l2)=str_L;
   S_out(w2+1:w2+r1,s2-r2+1:s2)=str_R;
   for ii=1:zn,
      S_out(w3+ii,l2+5+ii)='\'; 
   end;
   if ~isempty(str_crp),
      if (length(str_crp)+s2+8 < 80), 
         S_aux=repmat(' ',[s1,length(str_crp)+8]);
         S_aux(ceil(s1/2),9:end)=str_crp;
         S_out=[S_out,S_aux];
      else,
         S_out=char(S_out,['    ',str_crp]);
         if ls, S_out=char(S_out,' '); end; 
      end;
   end;
   S_out=char([iname,' = '],S_out);
   if ls, S_out=char(' ',S_out); end;
   disp(S_out);
else, 
   str_vr='';
   if ~isempty(F.frac.v), str_vr=[' in ',F.frac.v]; end;
   s1=F.frac.s(1); s2=F.frac.s(2);
   S_out=['  Left-denominator fraction (denominator\numerator)',str_vr,': ',int2str(s1),'-by-',int2str(s2)];
   if ~isempty(str_crp), S_out=char([S_out,','],['  ',str_crp]); end;
   if ls,
      S_out=char(' ',S_out,' ',[iname,'.denominator = '],str_L,...
         ' ',[iname,'.numerator = '],str_R); 
   else,
      S_out=char(S_out,[iname,'.denominator = '],str_L,...
         [iname,'.numerator = '],str_R);
   end;
   if strcmp(fr,'nice'), S_out=S_out(ls+2:end,:); end;
   disp(S_out);
end;

%end .. @ldf/display
