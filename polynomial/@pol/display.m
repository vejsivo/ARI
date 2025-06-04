function display(P)
%DISPLAY   Display polynomial matrix
%
% Syntax
%    DISPLAY(P)
% The display format is set by macro PFORMAT.

%       Author(s):  D. Henrion, S. Pejchova, M. Sebek 12-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 26-May-2000 10:21:27   $
%       &Revision: 3.0 &  $Date: 06-Sep-2000  S.Pejchova  $
%                         $Date: 21-Sep-2001  J.Jezek  

global PGLOBAL;
eval('PGLOBAL.FORMAT;', 'painit;');
fr = PGLOBAL.FORMAT;
ls = 1; orr = 0;

s1 = P.s(1); s2 = P.s(2);
vr = P.v; Pd = P.d; Pc = P.c;
if strcmp(vr,'z^-1'),
   vr1 = vr(1:3);
else, vr1 = [vr,'^'];
end;

if strcmp(get(0,'FormatSpacing'),'compact'),
   ls = 0;
end;
iname = [inputname(1),' = '];
if ls,
   iname = char('   ',iname,'   ');
end;
fr_b = get(0,'Format');
if strcmp(fr,'nice') & (size(Pc,3)<2),
   if ~size(Pc,3),
      Pc = zeros(s1,s2);
   end;
   eval([inputname(1),'=Pc(:,:)']);
elseif isempty(Pd),
   disp(iname);
   disp(['Empty polynomial matrix: ',int2str(s1),...
         '-by-',int2str(s2)]);
   if ls,
      disp(' '); 
   end;
elseif isinf(Pd),
   if ls,
      disp(' '); 
   end;
   disp(['Zero polynomial matrix: ',int2str(s1),...
         '-by-',int2str(s2),',  degree: -Inf']);
   disp([iname]);
   disp(zeros(s1,s2));
elseif Pd==0 & ~(strcmp(fr,'symbs')) & ~(strcmp(fr,'symb')), 
   if ls,
      disp(' ');
   end;
   disp(['Constant polynomial matrix: ',int2str(s1),...
         '-by-',int2str(s2)]);
   disp([iname]);
   if strcmp(fr,'symbr'),
      format rat;
   end;
   disp(Pc);
   set(0,'Format',fr_b);
else,
   switch fr,
   case {'coef','rcoef','block'},
      if ls,
         disp(' ');
      end;
      disp(['Polynomial matrix in ',vr,': ',int2str(s1),...
            '-by-',int2str(s2),',  degree: ',int2str(Pd)]);
      disp([iname]);
      if strcmp(PGLOBAL.ORDER,'reverse'),
         orr = 1;
      end;
      if orr & (strcmp(fr(end),'f')),
         fr = ['r',fr]; fr = strrep(fr,'rr','');
      end;
      
      switch fr,
      case 'coef',     
         for i = 0:Pd,
            disp(['  Matrix coefficient at ',vr1, int2str(i) ' :']);
            disp(Pc(:, :, 1+i));
         end;
      case 'rcoef',
         for i = Pd:-1:0,
            disp(['  Matrix coefficient at ',vr1, int2str(i) ' :']);
            disp(Pc(:, :, 1+i));
         end;
      case 'block',
         if orr,
            Pc = Pc(:,:,end:-1:1);
         end;
         disp(reshape(Pc, [s1 (Pd+1)*s2]));
      end; % switch fr
    
   case {'symb','symbs','symbr','rootr','rootc','nice'},
      switch fr,
      case 'symb',
         C = char(P);
      case {'symbs','nice'},
         C = char(P,2);  
      case 'symbr',
         C = char(P,'rat');
      case {'rootr','rootc'},
         C = char(P,fr);
      end;
      c1 = char(32*ones(s1,4));
      if strcmp(fr,'nice'),
         iname1=iname;
      elseif Pd==0,
         if ls,
            disp(' ');
         end;
         sTr_Au = ['Constant polynomial matrix: ',int2str(s1),...
            '-by-',int2str(s2)];
         iname1 = char([sTr_Au],iname); iname = iname1;
      else
         if ls,
            disp(' ');
         end;
         sTr_Au = ['Polynomial matrix in ',vr,': ',int2str(s1),...
            '-by-',int2str(s2),',  degree: ',int2str(Pd)];   
         iname1 = char([sTr_Au],iname);
      end;
      if (s1==1 & s2==1),
         disp(iname);
         disp([c1,C]); 
         if ls,
            disp(' ');
         end; 
      elseif s2==1,
         disp(iname);
         disp([c1,char(C{:,1})]);
         if ls,
            disp(' ');
         end;
      else,
         S = char(C{:,1}); S = [c1,S]; 
         col = 1; sumlg = size(S,2);
         for i = 2:s2+1,
            if i<=s2,
               Saux = char(C{:,i}); Saux = [c1,Saux];
               onelg = size(Saux,2);
            end;
            if (sumlg+onelg < 79) & (i<=s2),
               col = [col,i]; S = [S,Saux]; sumlg = sumlg+onelg;
            else,
               if length(col)~=s2,
                  if any(col==1),
                     disp(iname1);
                  end;
                  lgcol = length(col);
                  if lgcol==1,
                     disp(['  Column ',int2str(col)]);
                  else,
                     disp(['  Columns ',int2str(col(1)),...
                        ' through ',int2str(col(lgcol))]);
                  end;
               else,
                  disp(iname);
               end;
               disp(S);
               if ls,
                  disp(' ');
               end;
               col = i; S = Saux; sumlg = onelg;
            end;
         end; % for i=1:s2
      end; % if (s1==1...
   end; % switch rf
end;

%end .. @pol/display
