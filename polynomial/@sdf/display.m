function display(M)
%DISPLAY   Display scalar-denominator fraction
%
% Syntax
%    DISPLAY(M)
% The display format is set by macro PFORMAT.

%       Author(s):  S. Pejchova  15-02-00
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 16-Feb-2000 11:12:27   $
%                         $Date: 10-Apr-2000 11:01:00 S.Pejchova $
%                         $Date: 20-Apr-2000 18:00:00 J.Jezek  $
%                         $Date: 14-Sep-2000  S.Pejchova  $
%                         $Date: 14-Oct-2002  J.Jezek     $

global PGLOBAL
two=0;  ls=1;  str_crp=''; fr='';

eval('fr=PGLOBAL.FORMAT;',...
   'error(''Use PINIT to initialize the Polynomial Toolbox.'');');
if strcmp(get(0,'FormatSpacing'),'compact'), ls=0; end;
iname=inputname(1);
if strcmp(fr,'nice'),
   [aM1,cM1]=declass(M);
   if ~strcmp(cM1,'sdf'),
      eval([iname,'=aM1;  display(',iname,');']);
      return; 
   end; 
end;

if ~strcmp(M.frac.c,'cop?'),
   str_crp='coprime';
   if strcmp(M.frac.c,'ncop'),
      str_crp=['non',str_crp];
   end;
end;
if ~strcmp(M.frac.r,'red?'),
   s_r='reduced';
   if strcmp(M.frac.r,'nred'),
      s_r=['non',s_r];
   end;
   if isempty(str_crp), str_crp=s_r; 
   else, str_crp=[str_crp,', ',s_r]; 
   end;
end;
if ~strcmp(M.frac.p,'prop?'),
   s_r='proper';
   if strcmp(M.frac.p,'nprop'),
      s_r=['im',s_r];
   end;
   if isempty(str_crp), str_crp=s_r; 
   else, str_crp=[str_crp,', ',s_r]; 
   end;
end;
if strcmp(fr,'nice'), str_crp=''; end;

in_coef=strmatch(fr,['coef ';'rcoef';'block'],'exact');
if ~isempty(in_coef), two=1; end;

str_L=pdisp(M.frac.num);
fe1=findstr('Empt',str_L(ls+1,:));
if isempty(fe1), 
   fe1=findstr('matr',str_L(ls+1,:));
   if ~isempty(fe1), str_L(1:ls+1,:)=[]; end;
   if ~isempty(findstr('Column',str_L(ls+1,:))), two=1; end; 
else,
   str_L=[repmat(' ',[size(str_L,1),2]),str_L]; two=0;
end;   
str_L=str_L(:,1:max(find(any(str_L~=' ',1))));
if all(str_L(:,5)==' '), str_L=str_L(:,2:end); end;

str_R=pdisp(M.frac.den);
fe1=findstr('matr',str_R(ls+1,:));
if ~isempty(fe1), str_R(1:ls+1,:)=[]; end;
str_R=str_R(:,1:max(find(any(str_R~=' ',1))));
if all(str_R(:,5)==' '), str_R=str_R(:,2:end); end;

[l1,l2]=size(str_L);
[r1,r2]=size(str_R);
if (l1>20), two=1; end;

if ~two,
   s2=max(l2,r2);
   l_min=min(find(any(str_L~=' ',1)))-1;
   r_min=min(find(any(str_R~=' ',1)))-1;
   str_L=[repmat(' ',[l1,fix((s2-l2)/2)]),str_L];
   str_R=[repmat(' ',[r1,fix((s2-r2)/2)]),str_R];
   str_zn=['    ',repmat('-',[1,max(0,s2-4)])];
   S_out=strvcat(str_L,str_zn,str_R); 
   s1=l1+r1+1; 
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
   if ~isempty(M.frac.v), str_vr=[' in ',M.frac.v]; end;
   s1=M.frac.s(1); s2=M.frac.s(2);
   S_out=['  Scalar-denominator fraction',str_vr,': ',int2str(s1),'-by-',int2str(s2)];
   if ~isempty(str_crp), S_out=char([S_out,','],['  ',str_crp]); end;
   if ls,
      S_out=char(' ',S_out,' ',[iname,'.numerator = '],str_L,...
         ' ',[iname,'.denominator = '],str_R); 
   else,
      S_out=char(S_out,[iname,'.numerator = '],str_L,...
         [iname,'.denominator = '],str_R); 
   end;
   if strcmp(fr,'nice'), S_out=S_out(ls+2:end,:); end;
   disp(S_out);
end;

%end .. @sdf/display
