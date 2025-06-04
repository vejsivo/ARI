function Str = pdisp(P)
%PDISP  Display of a polynomial matrix without printing the name.
%
% Syntax
%    PDISP(P)
% it's the same as DISPLAY(P) or as leaving the semicolon off an
% expression without printing the matrix name.  
% The display format is set by macro PFORMAT.
% If output argument is present
%    STR = PDISP(P)
% then the displaying string is saved in STR for later using
% and DISP(STR) is the same as PDISP(P).

%       Author(s):  S. Pejchova 16-07-99
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 29-May-2000 15:20:27   $
%                         $Date: 12-Sep-2000  S.Pejchova  $
%                         $Date: 31-Oct-2000  S.Pejchova  $
%                         $Date: 07-Jan-2002  J.Jezek     $

global PGLOBAL;
eval('PGLOBAL.VARIABLE;', 'painit;');

ls=1; S_out='';
fr=PGLOBAL.FORMAT;
narginchk(1,1);
% error(nargchk(1,1,nargin));	%REMOVED IN NEW MATLABS
eval('P=pol(P);', 'error(peel(lasterr));');
vr=P.v; Pd=P.d; [s1,s2]=size(P); Pc=P.c;
if strcmp(vr,'z^-1'),
   vr1 = vr(1:3);
else, vr1 = [vr,'^'];
end;
if strcmp(get(0,'FormatSpacing'),'compact'), ls=0; end;

fr_b=get(0,'Format');
if ~Pd & (strcmp(fr,'coef')|strcmp(fr,'rcoef')|strcmp(fr,'block')),
   fr='symb';
end;

if isempty(Pd),
   S_out=['  Empty'];
   if ~strcmp(fr,'nice'),  S_out=[strrep(S_out,' ',''),' polynomial']; end;
   S_out=[S_out,' matrix: ',int2str(s1),'-by-',int2str(s2)];
   if ls, S_out=char(' ',S_out,' '); end;
elseif isinf(Pd),
   x_str=['Zero polynomial matrix: ',int2str(s1),'-by-',int2str(s2),...
         ',  degree: -Inf'];
   S_out=repmat('     0',[s1, s2]);
   if ls, S_out=char(' ',S_out,' '); x_str=char(' ',x_str); end;
   if ~strcmp(fr,'nice'), S_out=char(x_str,S_out); end;
else,
   n_lg=[];
   switch fr,
   case {'symbs','nice'}, n_lg=2;  
   case 'symbr', n_lg='rat';
   case {'rootr','rootc'}, n_lg=fr;
   end;
   switch fr,
   case {'coef','rcoef','block'},
      S_out=['Polynomial matrix in ',vr,': ',int2str(s1),'-by-',int2str(s2),...
             ',  degree: ',int2str(Pd)];
      if ls, S_out=char(' ',S_out); end;
      if strcmp(fr,'block'),
         Pc=Pc(:,:);
         Sii=m2str(Pc,n_lg);
         if ls, 
            S_out=char(S_out,' ',Sii); 
         else,
            S_out=char(S_out,Sii);
         end; 
      else, 
         jj=0:Pd;
         if strcmp(fr,'rcoef'), jj=Pd:-1:0; end;
         for ii=jj,
            Sii=m2str(Pc(:, :, 1+ii),n_lg);
            Str_ii=['  Matrix coefficient at ',vr1, int2str(ii) ' :'];
            if ls, 
               S_out=char(S_out,' ',Str_ii,' ',Sii); 
            else,
               S_out=char(S_out,Str_ii,Sii);
            end; 
         end;
      end; 
     
   case {'symb','symbs','symbr','rootc','rootr','nice'},
      if isempty(n_lg), C=char(P); 
      else, C=char(P,n_lg);
      end;
      c1=repmat([' '],[s1,4]); S_out='';
      if Pd<0,
         iname1=['Zero polynomial matrix: ',int2str(s1),'-by-',int2str(s2),...
               ',  degree: -Inf'];
      elseif ~Pd,
         iname1=['Constant polynomial matrix: ',int2str(s1),'-by-',int2str(s2)];
      else,
         iname1=['Polynomial matrix in ',vr,': ',int2str(s1),...
                 '-by-',int2str(s2),',  degree: ',int2str(Pd)];   
      end;
      if ls,  iname2=char(' ',iname1); iname1=char(' ',iname1,' '); 
      else, iname2=iname1; 
      end;
      if s2==1,
         if s1==1, S_out=[c1,C]; 
         else,  S_out=[c1,char(C{:,1})]; 
         end;
         if ls, S_out=char(' ',S_out,' '); end;
         if ~Pd&(~strcmp(fr,'nice')), S_out=char(iname2,S_out); end;
      else,
         S=char(C{:,1}); S=[c1,S]; 
         col=1; sumlg=size(S,2);
         for i=2:s2+1,
            if i<=s2,
               Saux=char(C{:,i}); Saux=[c1,Saux];
               onelg=size(Saux,2);
            end;
            if (sumlg+onelg < 79)&(i<=s2),
               col=[col,i]; S=[S,Saux]; sumlg=sumlg+onelg;
            else,
               if length(col)~=s2,
                  if any(col==1)&(~strcmp(fr,'nice')), S_out=iname1; end;
                  lgcol=length(col);
                  S_col=['  Column ',int2str(col(1))]; 
                  if lgcol~=1,
                     S_col=strrep(S_col,'n','ns');
                     S_col=[S_col,' through ',int2str(col(lgcol))];
                  end;
                  if isempty(S_out), S_out=S_col; 
                  else,  S_out=char(S_out,S_col); 
                  end;
               elseif ~Pd&(~strcmp(fr,'nice')), S_out=iname1; end;
               if isempty(S_out)&~ls, 
                  S_out=S;
               else,
                  S_out=char(S_out,S); 
               end;
               if ls, S_out=char(S_out,' '); end;
               col=i; S=Saux; sumlg=onelg;
            end;
         end; % for i=1:s2
      end; % if (s2==1...
   end; % switch rf
end;
if strcmp(fr,'nice')&ls&any(S_out(1,:)~=' '),  S_out=char(' ',S_out); end;
if ~nargout,
   disp(S_out);
else,
   Str=S_out; 
end;
   
% SUBFUNCTION - converts the Matlab matrix A into string representation
%               with respect to defined format or number of precision (fmt),

function sr=m2str(A,fmt),
[v1,v2]=size(A);
if isempty(fmt), chA=char(pol(A));
else, chA=char(pol(A),fmt);
end;
if v1*v2==1,                             % added by J.Jezek 07-Jan-2002
   schA=chA; chA=cell(1,1); chA{1}=schA; %
end;                                     %
S1=[repmat(' ',[v1*v2,4]),char(chA')]';
sr=reshape(S1(:),[prod(size(S1))/v1,v1]);
sr=sr';

%end .. m2str

%end .. pdisp
