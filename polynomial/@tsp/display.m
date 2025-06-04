function display(T)
%DISPLAY   Display two-sided polynomial matrix
%
% Syntax
%    DISPLAY(T)
% The display format is set by macro PFORMAT.

%       Author(s):  S. Pejchova  08-7-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 17-Sep-1999 11:52:27   $
%       $Revision: 3.0 $  $Date: 30-May-2000 15:14:27   $
%       &Revision: 3.0 &  $Date: 07-Sep-2000  S.Pejchova  $
%                         $Date: 24-Sep-2001  J.Jezek   $

global PGLOBAL;
eval('PGLOBAL.FORMAT;', 'painit;');
fr = PGLOBAL.FORMAT;
ls = 1; orr = 0;

s1 = T.s(1); s2 = T.s(2);
Tp = T.p; Tc = Tp.c; offst = T.o; Tr = T.r;
vr = 'z^-1'; vr1 = 'z^-';

if strcmp(get(0,'FormatSpacing'),'compact'),
   ls = 0;
end;
iname = [inputname(1),' = '];
if ls,
   iname = char('   ',iname,'   ');
end;
fr_b = get(0,'Format');
if strcmp(fr,'nice') & (size(Tc,3)<2) & (~offst),
   if ~size(Tc,3),
      Tc = zeros(s1,s2);
   end;
   eval([inputname(1),'=Tc(:,:)']);
elseif isempty(Tr),
   disp(iname);
   disp(['Empty two-sided polynomial matrix: ',int2str(s1),...
         '-by-',int2str(s2)]);
   if ls,
      disp(' ');
   end;
elseif isinf(Tr),
   if ls,
      disp(' ');
   end;
   disp(['Zero two-sided polynomial matrix: ',int2str(s1),...
         '-by-',int2str(s2),',  degree: -Inf,  tdegree: +Inf']);
   disp([iname]);
   disp(zeros(s1,s2));
elseif (~Tr) & (~offst) & ~(strcmp(fr,'symbs')) & ~(strcmp(fr,'symb')), 
   if ls,
      disp(' ');
   end;
   disp(['Constant two-sided polynomial matrix: ',int2str(s1),...
      '-by-',int2str(s2)]);
   disp([iname]);
   if strcmp(fr,'symbr'),
      format rat;
   end;
   disp(Tp{0});
   set(0,'Format',fr_b);
else,
   Tr_str = [int2str(Tr+offst),',  tdegree: ',int2str(offst)];
   
   switch fr,
   case {'coef','rcoef','block'},
      if ls,
         disp(' ');
      end;    
      disp(['Two-sided polynomial matrix in z: ',int2str(s1),...
            '-by-',int2str(s2),',  degree: ',Tr_str]);
      disp([iname]);
      if strcmp(PGLOBAL.ORDER,'reverse'),
         orr = 1;
      end;
      if orr & (strcmp(fr(end),'f')),
         fr = ['r',fr]; fr = strrep(fr,'rr','');
      end;
         
      switch fr,
      case 'coef',
         for i = 1:(Tr+1),
            disp(['  Matrix coefficient at z^', int2str(i+offst-1) ' :']);
            disp(Tc(:, :, i));
         end;          
      case 'rcoef',
         for i = (Tr+1):-1:1,
             disp(['  Matrix coefficient at z^', int2str(i+offst-1) ' :']);
             disp(Tc(:, :, i));
         end;             
      case 'block',
         if orr,
            Tc = Tc(:,:,end:-1:1);
         end;
         disp(Tc(:,:));
      end; % switch fr
    
   case {'symb','symbs','symbr','rootr','rootc','nice'},
      switch fr,
      case 'symb', 
         C = char(T);
      case {'symbs','nice'},
         C = char(T,2);  
      case 'symbr',
         C = char(T,'rat');
      case {'rootr','rootc'},
         C = char(T,fr);          
      end;
      c1 = char(32*ones(s1,4));
      if strcmp(fr,'nice'),
         iname1 = iname;
      elseif ~Tr & ~offst,
         if ls,
            disp(' ');
         end;
         sTr_Au = ['Constant two-sided polynomial matrix: ',int2str(s1),...
               '-by-',int2str(s2)];
         iname1 = char([sTr_Au],iname); iname = iname1;
      else,
         if ls,
            disp(' ');
         end;
         sTr_Au = ['Two-sided polynomial matrix in z: ',int2str(s1),...
                 '-by-',int2str(s2),',  degree: ',Tr_str];   
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

%end .. @tsp/display
