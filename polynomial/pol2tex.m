function tex_str=pol2tex(varargin)
%POL2TEX  LaTeX representation of a polynomial object.
%
% The commmand
%    TEX_STR =  POL2TEX(A1,A2,..,AN) 
% where A1,..,AN are standard MATLAB constant matrices or polynomial matrices,
% returns a string TEX_STR with a LaTeX control sequences to produce these matrices
% in a array environment of the LaTeX mathematical mode.
% 
%    TEX_STR =  POL2TEX(A1,A2,..,AN,'FILE_NAME')
% appends TEX_STR with LaTeX representation of matrices A1,..,AN to existing TEX file
% FILE_NAME.TEX. If the FILE_NAME.TEX does not exist, then TEX_STR is saved in a
% new-created TEX file.

%       Author:  S. Pejchova 17-6-99
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 01-Jun-2000 11:34:27   $
%                         $Date: 10-May-2003  S.Pejchova  $

global PGLOBAL;
eval('PGLOBAL.FORMAT;', 'painit;');
ls=1; file_nm=[]; tex_str='';
fr=PGLOBAL.FORMAT;

ni=nargin;
if ni==0,
   error('Not enough input arguments.');
end;
fr_b=get(0,'Format');
for ii=1:ni,
   arg_x=varargin{ii};
   if ischar(arg_x),
      tx1=findstr(arg_x,'.tex');
      if ~isempty(tx1)
         file_nm=arg_x;
      else,
         file_nm=[arg_x,'.tex'];
      end;
   else,
      tex_str1='';
      if ~isempty(inputname(ii)), tex_str1=[sprintf('\n'),inputname(ii),'=']; end;
      cl_argx=class(arg_x);
      Ind_cl=strmatch(cl_argx,{'pol';'double';'tsp'},'exact');
   if isempty(Ind_cl),
         % Matrix is not double or polynomial
         error('Invalid class of input argument.');
   elseif isempty(arg_x),
         % Empty double or polynomial matrices  
         arg_x=arg_x(:,:); 
         ep1=max(size(arg_x),[1,1]);
         if any(ep1~=1),
            t_sz=mx2tex(ones(ep1));
            t_sz=strrep(t_sz,'1',' \; ');
         else,
            t_sz='[ \; ]';
         end;
         tex_str=[tex_str,sprintf('\n'),'$$',tex_str1,sprintf('\n'),...
                  t_sz,sprintf('\n'),'$$'];
         
   elseif strcmp(cl_argx,'double'),
         % Standard Matlab double matrices
         tex_str=[tex_str,sprintf('\n'),'$$',tex_str1,sprintf('\n'),...
                  mx2tex(arg_x),sprintf('\n'),'$$'];
   else,
      [s1,s2]=size(arg_x); Po=0;
      if strcmp(cl_argx,'pol'), Pd=arg_x.d; Pm=Pd;
      else, Pd=arg_x.r; Po=arg_x.o; Pm=arg_x.d;
      end;
      
         if (Pd<=0)&(Po==0),
            % Zero and constant polynomial matrices
            tex_str=[tex_str,sprintf('\n'),'$$',tex_str1,sprintf('\n'),...
                  mx2tex(arg_x{0},fr),sprintf('\n'),'$$'];
         else, 
            vr=strrep(arg_x.v,'-1','{-1}');  vr1=vr(1);
            switch fr,
            case 'block',
             tex_str=[tex_str,sprintf('\n'),'$$',tex_str1,sprintf('\n'),...
                  mx2tex(arg_x{:}),sprintf('\n'),'$$'];
              
             case {'coef','rcoef'},
               t_sz='';  t_ind=Po:Pm;
               if strcmp(fr,'rcoef'), t_ind=fliplr(t_ind); end;
               for kk=t_ind,
                  nk=kk;
                  if length(vr)>1, nk=-kk; end;
                  t_sz1=[sprintf('\n'),vr1,'^{',int2str(nk),'}'];
                  if ~kk, t_sz1='';, elseif kk==1, t_sz1=[sprintf('\n'),vr];, end
                  if isempty(t_sz), 
                     t_sz=[mx2tex(arg_x{kk}),t_sz1];
                  else,
                     t_sz=[t_sz,sprintf('\n'),'+',sprintf('\n'),mx2tex(arg_x{kk}),t_sz1];
                  end;
               end;
               tex_str=[tex_str,sprintf('\n'),'$$',tex_str1,sprintf('\n'),...
                        t_sz,sprintf('\n'),'$$'];
             otherwise,
               switch fr,
               case 'symb', C=char(arg_x);
               case {'symbs','nice'}, C=char(arg_x,2);
               case 'symbr',C=char(arg_x,'rat');
               case {'rootr','rootc'}, C=char(arg_x,fr);
               end;
               ent_str='';
               for ix=1:s1,
                  for jx=1:s2,
                    if all([s1,s2]==1), xx_str=C;
                    else, xx_str=C{ix,jx};
                    end;
                    xx_str=strrep(xx_str,' ','');
                    f_e=fliplr(sort([findstr(xx_str,'e-'),findstr(xx_str,'e+')])); 
                    for k=f_e,
                       rn_str=['*10^{',int2str(eval(xx_str(k+1:k+4))),'} '];
                       xx_str=[xx_str(1:k-1),rn_str,xx_str(k+5:end)];
                    end;
                    vr2=[vr1,'^']; f_v=[];
                    if length(vr2)<length(xx_str),f_v=fliplr(findstr(xx_str,vr2)); end;
                    for k=f_v,
                       v_str=xx_str(k+3:end);
                       ff_v=sort([findstr(v_str,'+'),findstr(v_str,'-'),...
                             findstr(v_str,'('),findstr(v_str,')')]);
                       if isempty(ff_v),
                          xx_str=[xx_str(1:k+1),'{',xx_str(k+2:end),'}'];
                       else,
                          xx_str=[xx_str(1:k+1),'{',xx_str(k+2:(k+1+ff_v(1))),...
                                '}',xx_str((k+2+ff_v(1)):end)];
                       end;
                    end;
                    if strcmp(fr,'symbr')|(strcmp(fr,'symb')&strcmp(fr_b,'rational'))|...
                          (strncmp(fr,'root',4)&strcmp(fr_b,'rational')),

                       ex_s1=cumsum((xx_str=='(')-(xx_str==')'));
                       ex_s2=cumsum((xx_str=='{')-(xx_str=='}'));
                       if strncmp(fr,'root',4),
                          ex_s3=((xx_str=='-')|(xx_str=='+')) &...
                             (ex_s1<2) & (ex_s2==0);
                       else,
                          ex_s3=((xx_str=='-') | (xx_str=='+')) & (ex_s1==0) & (ex_s2==0);
                       end;
                       %ex_s3=((xx_str=='-') | (xx_str=='+')) & (ex_s1==0) & (ex_s2==0);
                       fe2=find(ex_s3(2:end));
                       fe2=[0,fe2]; ex_s0=xx_str; xx_str='';
                       for nn=fliplr(fe2),
                          ex_s1=ex_s0(nn+1:end);
                          f_0=sort([findstr(ex_s1,vr1),findstr(ex_s1,['*',vr1])]); 
                          ex_0end=''; ex_0bg=ex_s1;
                          if ~isempty(f_0), 
                             ex_0bg=ex_s1(1:(f_0(1)-1));
                             ex_0end=ex_s1(f_0(1):end);
                             if strncmp(ex_0end,'*',1), ex_0end=ex_0end(2:end); end;
                          end;
                          f_1=find(ex_0bg=='/');
                          if ~isempty(f_1), 
                             ex_0bg=strrep(ex_0bg,'/','}{');
                             if strncmp(ex_0bg,'-',1),
                                ex_0bg=['- \frac{',ex_0bg(2:end),'}'];
                             elseif strncmp(ex_0bg,'+',1),
                                ex_0bg=['+ \frac{',ex_0bg(2:end),'}'];
                             else,
                                ex_0bg=[' \frac{',ex_0bg,'}'];
                             end;
                             if strcmp(ex_0bg(end-1:end),')}'), ex_0bg(end-1:end)='})'; 
                             elseif strcmp(ex_0bg(end-2:end),')(}'), ex_0bg(end-2:end)='})(';
                             end
                          end;
                          ex_s0(nn+1:end)=[];
                          xx_str=[ex_0bg,ex_0end,xx_str];
                       end;
                    end;        
                    
                    if ~strncmp(xx_str,'-',1), xx_str=['\;\;\;',xx_str]; end;
                    if jx>1, xx_str=[' & ',xx_str];  end
                    ent_str=[ent_str,xx_str];
                 end; % for jx
                 if ix~=s1, 
                    ent_str=[ent_str,' \\'];
                    if strcmp(fr,'symbr'),ent_str=[ent_str,' \\'];, end
                 end;
                 ent_str=[ent_str,sprintf('\n')];
               end; % for ii
               if any([s1,s2]>1),
                  ent_str=['\left[ \begin{array}{',char(108*ones(1,s2)),'}',...
                           sprintf('\n'),ent_str];
                  ent_str=[ent_str,'\end{array} \right]'];
               end;       
               tex_str=[tex_str,sprintf('\n'),'$$',tex_str1,sprintf('\n'),...
                        ent_str,sprintf('\n'),'$$'];
              
            end; %switch fr
            
         end; % any(isinf(Pd))
         
      end; % isempty(Ind_cl)
      
          
         
   end; % if ischar(arg_x),
end; % for ii=1:ni,
if ~isempty(tex_str), tex_str=[tex_str,sprintf('\n')]; end;

if ~isempty(file_nm)&~isempty(tex_str),
   f_id1=fopen(file_nm,'at');
   fprintf(f_id1,'%s',tex_str);
   fclose(f_id1);
end;



%------------SUBFUNCTION---------------------------
function out_str=mx2tex(A,form_A),
%LaTeX representation of nonempty standard Matlab 2-D matrix
%without the name (similar as DISP function)
fA=1;   [s1,s2]=size(A);
if nargin==1,  form_A=get(0,'Format'); end;
if strcmp(form_A,'rational')|strcmp(form_A,'symbr'), 
   out_str='';
elseif strncmp(form_A,'long',4), 
   out_str=mat2str(A,15);
elseif strcmp(form_A,'symbs')|strcmp(form_A,'bank'),
   out_str=mat2str(A,2);
else,
   out_str=mat2str(A,4);
end;
if isempty(out_str),
   for ii=1:s1,
      for jj=1:s2,
         [n1,d1]=rat(A(ii,jj));
         sxN = num2str(abs(n1));
         if abs(d1)~=1, sxN = ['\frac{',sxN,'}{',int2str(d1),'}']; end;
         if n1>0,
            sxN=[' \;\;\; ',sxN];
         else,
            sxN=['-',sxN];
         end;
         if jj>1, sxN=[' & ',sxN];  end
         out_str=[out_str,sxN];
      end; % for jj
      if ii~=s1, out_str=[out_str,' \\ \\'];, end;
      out_str=[out_str,sprintf('\n')];
   end; % for ii
else,
   out_str=strrep(out_str,' -','&-');
   out_str=strrep(out_str,' ','&?');
   out_str=strrep(out_str,';-',[' \\',sprintf('\n'),'-']);
   out_str=strrep(out_str,';',[' \\',sprintf('\n'),'?']);
   out_str=strrep(out_str,'[-','-');
   out_str=strrep(out_str,'[','?');
   out_str=strrep(out_str,'?',' \;\;\; ');
   out_str=strrep(out_str,'&',' & ');
   out_str=strrep(out_str,']',sprintf('\n'));
   f_e=sort([findstr(out_str,'e-'),findstr(out_str,'e+')]); 
   if ~isempty(f_e),
      f_e=fliplr(f_e);
      for k=f_e,
         rn_str=['*10^{',int2str(eval(out_str(k+1:k+4))),'} '];
         out_str=[out_str(1:k-1),rn_str,out_str(k+5:end)];
      end;
   end;
end;
if length(A)>1,
   out_str=['\left[ \begin{array}{',char(108*ones(1,s2)),'}',sprintf('\n'),out_str];
   out_str=[out_str,'\end{array} \right]'];
end;

%end .. mxtex

%end .. pol2tex
