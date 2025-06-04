function pme(command_str)
%PME  Polynomial matrix editor
%
% The Polynomial Matrix Editor is a GUI designed to create and modify 
% polynomial matrices (POL objects) easily and interactively. It handles 
% standard Matlab matrices (2-D DOUBLE objects) as well.
%
% QUICK START:
%
% Type PME to open a Polynomial Matrix Editor window.
% This window displays existing polynomial matrices and allows you to
% create a new one.
%
% * To create a new polynomial matrix, type its name and size in the
%   first (editable) line and then click 'Open'.
%
% * To modify an existing polynomial matrix while keeping its size and
%   other properties, just find it in the list and then double click 
%   on the corresponding row.
%
% * To modify an existing polynomial matrix substantially (changing also
%   its name, size, variable symbol, etc.), first click on the 
%   corresponding row in the list to move it up, then type in the new 
%   required properties and finally click 'Open'.
%
% Each of the actions above opens a Matrix Pad to edit the matrix. 
%
% On the Matrix Pad, the matrix entries apper in boxes. If an entry is 
% too long and cannot be fully displayed, its box takes a slightly different 
% color (usually more pinkish). 
%    
% To edit an entry, just click on its box. Then its box becomes editable 
% and large enough to see the whole content. You can write in anything that 
% follows MATLAB syntax and results in a polynomial matrix or scalar (so 
% you can use existing MATLAB variables, functions, etc.). The program is 
% even more intelligent and allows some notation going beyond the MATLAB 
% syntax. For instance, you can drop the * (times) operator between a 
% coefficient and the related polynomial symbol (provided that the coefficient
% comes first) so you can type 2s as well as 2*s. You can also type zi as 
% well as z^-1. To complete editing of the entry, press the 'Entry' key or 
% close the box by mouse. If you have entered an expression that cannot be 
% processed then an error is reported and the original box content is recovered.
%
% To bring the new or modified matrix into the MATLAB workspace click 'Save' or
% 'Save As'.
%
% MORE DETAILS:
%
% POLYNOMIAL MATRIX EDITOR WINDOW - BUTTONS
%
% * 'Open' Button: Clicking the 'Open' button opens a new Matrix Pad for the 
%   matrix specified by the first (editable) line. At most 4 Matrix Pads can be 
%   open at the same time.
%
% * 'Refresh' Button: Clicking the 'Refresh' button updates the list of the main 
%   window to reflect possible changes in the workspace since opening the PME 
%   or the last time it was refreshed.
%
% * 'Close' Button: Clicking the 'Close' button terminates the current PME session.
%
% POLYNOMIAL MATRIX EDITOR WINDOW - OPTIONS
%    
% * Menu 'List': By using the menu 'List' you can set the class of objects 
%   displayed by the main PME window. It may be POL objects (default) or 2-D 
%   DOUBLE objects (standard MATLAB matrices), or both POL and DOUBLE matrices 
%   at the same time.
%
% * Menu 'Symbol': By using the menu 'Symbol' you can choose the symbol (s,z,...) 
%   that is automatically offered for the 'New' matrix item of the list.
%
% MATRIX PAD - BUTTONS
%
% * 'Save' Button: Clicking 'Save' brings the Matrix Pad contents into the MATLAB 
%   workspace.
%
% * 'Save As' Button: Clicking 'Save As' puts the matrix on the Pad into the MATLAB 
%   workspace under another name. 
%
% * 'Browse' Button: Clicking 'Browse' moves the cursor to the main PME window (the 
%    same as directly clicking there).
%
% * 'Close' Button: Clicking 'Close' closes the Matrix Pad window.

%       Author(s):  S. Pejchova, M. Sebek 28-8-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 27-Apr-1999 17:42:34   $

global PGLOBAL;
eval('PGLOBAL.VARIABLE;', 'painit;');
vvx = PGLOBAL.VARIABLE;

if exist('pmexx2.mat')==2,
   delete pmexx2.mat;
end
if nargin==0
   command_str='initialize';
end
if ~ischar(command_str)|isempty(command_str), error('Illegal input argument.'); end;

switch command_str
   %                                                     INITIALIZATION
   case 'initialize',
     s_figs=get(0,'children');
     fig_ex=0;
     for fgr=s_figs'
       fig_ex=strcmp(get(fgr,'name'),'Polynomial Matrix Editor');
       if fig_ex
          figure(fgr);, return;
       end
     end
     pos_struct=struct('strlg',0,'s1',0,'s2',0,'s3',0,'char1',1,'ls',0,'zd',0');
     m_fig=figure('numbertitle','off','name','Polynomial Matrix Editor',...
                  'units','normalized','MenuBar','none',...
                  'ResizeFcn','pme(''cb_rsf2'');',...  
                  'UserData',pos_struct,...
                  'position',[0.4 0.53 0.55 0.37],'Tag','Small_PmE');

     %                                            Figure Frame
     uicontrol(m_fig,'style','frame','units','normalized',...
                   'BackgroundColor',[0.5 0.55 0.7],...
                   'ForegroundColor',[1 1 1],'Tag','S_frfig1',...
                   'position',[0.01 0.01 0.98 0.98]);
     %                                            Menu-Class Browser
     inp_fl=[];          
     menu1=uimenu(m_fig,'Label','List','Tag','S_pop1');
     uimenu(menu1,'Label','POL','Checked','on','CallBack',...
                 ['if exist(''pmexx2.mat'')==2,','load pmexx2; end;',...
                  'save pmexx1;','pme(''cb_type1'')']);
     uimenu(menu1,'Label','DOUBLE','Checked','off','CallBack',...
                 ['if exist(''pmexx2.mat'')==2,','load pmexx2; end;',...
                  'save pmexx1;','pme(''cb_type2'')']);
     uimenu(menu1,'Label','POL && DOUBLE','Checked','off','CallBack',...
                 ['if exist(''pmexx2.mat'')==2,','load pmexx2; end;',...
                    'save pmexx1;','pme(''cb_type3'')']);
              
     %                                             Menu-Global Symbol
     val_var_str=['s','d','p','q','z','z^-1'];         
     val_var=findstr(vvx,val_var_str); val_var=val_var(1);
     menu2=uimenu(m_fig,'Label','Symbol','UserData',vvx,'Tag','S_var1');
     for ii=1:6,
        if ii==val_var, sii='on'; else, sii='off'; end
        if ii==6, s_var='z^-1'; else, s_var=val_var_str(ii); end 
        uimenu(menu2,'Label',s_var,'Checked',sii,'CallBack','pme(''cb_var'')');
     end;   
        
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
     s_l1=uicontrol('Parent',m_fig,'Style','listbox','units','points',...
                   'Position',[20 31 290 64.5],...
                   'FontName','Courier','FontSize',9,...
                   'Visible','off','Enable','Off','string','N');
     xn1=get(s_l1,'Extent');
     set(s_l1,'String','Ne');    xn2=get(s_l1,'Extent');
     pos_l1=get(s_l1,'Position');
     delete(s_l1);
     char1=xn2(3)-xn1(3); ls=ceil(xn1(4));
     zd=xn1(3)-char1; 
     strlg=fix((pos_l1(3)-zd-ls)/char1);
     ps3=strlg*char1+zd+ls+1;
     pos_l1(3)=ps3;
     ps1=pos_l1(1); ps2=pos_l1(2); ps4=pos_l1(4);
     s3=strlg-3; s2=strlg-9-fix((2*zd+2+ls)/char1);
     s1=strlg-22-fix((3*zd+4+ls)/char1);
     pos_struct.strlg=strlg; pos_struct.s1=s1;
     pos_struct.s2=s2;       pos_struct.s3=s3;
     pos_struct.char1=char1; pos_struct.zd=zd;
     pos_struct.ls=ls;
     set(m_fig,'UserData',pos_struct);
     
     va4=round(ls); 
     va3=4*char1+ls+zd; va2=ps2+ps4+2; va1=ps1+ps3-1-va3;
           
     ca3=6*char1+ls+zd; ca1=ps1+ps3-3-va3-ca3;
     
     sz3=13*char1+zd;   sz1=ps1+ps3-5-va3-ca3-sz3;
     
     na3=ps3-va3-ca3-sz3-7;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %                                            Text-Disp.Class
     uicontrol('Parent',m_fig,'style','text','units','points',...
                   'BackgroundColor',[0.5 0.55 0.7],...
                   'ForegroundColor',[1 1 1],...
                   'FontSize',9,'FontWeight','demi', ...
                   'HorizontalAlignment','left',...
                   'Position',[20 136.5 ps3 ls], ...
                   'String','', ...
                   'Tag','S_type1',...
                   'UserData',inp_fl);
                

     %                                            Frame-table
     uicontrol('Parent',m_fig,'Style','frame','Units','points',...
                   'BackgroundColor',[1 1 1],...
                   'ForegroundColor',[1 1 1],'Tag','S_frta1',...
                   'Position',(pos_l1+[-2, -2, 4, 8+2*va4]));
                
     %                                            Text-Name
     uicontrol('Parent',m_fig,'style','text','units','points',...
                   'Position',[ps1, va2+va4+2, na3, va4], ...
                   'BackgroundColor',[0.5 0.55 0.7],...
                   'ForegroundColor',[1 1 1],...
                   'FontSize',9,'FontWeight','demi', ...
                   'HorizontalAlignment','center',...
                   'Tag','S_name1','String','Name');
     %                                            Text-Size
     uicontrol('Parent',m_fig,'Style','text','Units','points', ...
	                'Position',[sz1, va2+va4+2, sz3 ,va4], ...
                   'BackgroundColor',[0.5 0.55 0.7], ...
                   'ForegroundColor',[1 1 1],...
	                'FontSize',9,'FontWeight','demi', ...
	 	             'String','Size','UserData',[],...
	                'Tag','S_size1');
     %                                            Text-Class
     uicontrol('Parent',m_fig,'Style','text','Units','points', ...
	                'Position',[ca1, va2+va4+2, ca3, va4], ...
                   'BackgroundColor',[0.5 0.55 0.7], ...
	                'FontSize',9,'FontWeight','demi', ...
	                'ForegroundColor',[1 1 1], ...
	 	             'String','Object', ...
	                'Tag','S_class1');
     %                                            Text-Variable
     uicontrol('Parent',m_fig,'Style','text','Units','points', ...
	                'Position',[va1, va2+va4+2, va3, va4], ...
                   'BackgroundColor',[0.5 0.55 0.7], ...
	                'FontSize',9,'FontWeight','demi', ...
	                'ForegroundColor',[1 1 1], ...
	 	             'String','Variable', ...
                   'Tag','S_symb2');
               
     %                                            Pop-Class
     uicontrol('Parent',m_fig,'Style','popupmenu','Units','points',...
                   'Position',[ca1, va2, ca3, va4],...
                   'BackgroundColor',[0.9 0.9 1],...
                   'ForegroundColor',[0 0 0],...
                   'HorizontalAlignment','left', ...
                   'FontName','Courier','FontSize',9,...
                   'Tag','S_pop3','Enable','On',...
                   'String',['pol   ';'double'],'Value',1,...
                   'CallBack','pme(''cb_pop'')');
     %                                            Pop-Variable
     indeter_str=char('s','d','p','q','z','z^-1','');
     uicontrol('Parent',m_fig,'Style','popupmenu','Units','points',...
                   'Position',[va1, va2, va3, va4],...
                   'BackgroundColor',[0.9 0.9 1],...
                   'ForegroundColor',[0 0 0],...
                   'HorizontalAlignment','left', ...
                   'FontName','Courier','FontSize',9,...
                   'Tag','S_pop2','Enable','On',...
                   'String',indeter_str,...
                   'Value',val_var);
                
     %                                            Edit-Name
     uicontrol('Parent',m_fig,'Style','edit','Units','points',...
                   'Position',[ps1, va2-1, na3,va4+1],...
                   'BackgroundColor',[0.9 0.9 1],...
                   'ForegroundColor',[0 0 0],...
                   'HorizontalAlignment','left', ...                
                   'FontName','Courier','FontSize',9,...
                   'Tag','S_edit1','Enable','On',...
                   'UserData',[],'string','New');
                 
                
     %                                            Edit-Size
     uicontrol('Parent',m_fig,'Style','edit','Units','points',...
                   'Position',[sz1, va2-1, sz3, va4+1],...
                   'BackgroundColor',[0.9 0.9 1],...
                   'ForegroundColor',[0 0 0],...
                   'HorizontalAlignment','left', ...
                   'FontName','Courier','FontSize',9,...
                   'Tag','S_edit2','Enable','On',...
                   'string','3-by-3','UserData',[]);
                
                   
     %                                            List...
     StR0=char(32*ones(1,strlg));  val_list=1;
     user_struct=struct('name','New','size','3-by-3','type','pol','symb','','matr',pol([]));
     uicontrol('Parent',m_fig,'Style','listbox','units','points',...
                   'Position',[ps1, ps2, ps3, ps4],...
                   'BackgroundColor',[1 1 1],...
                   'ForegroundColor',[0 0 0],...
                   'FontName','Courier','FontSize',9,...
                   'Tag','S_list1','Enable','Off',...
                   'string',StR0,...
                   'value',val_list,...
                   'UserData',user_struct,...
                   'CallBack','pme(''cb_list'')');
                
                
     %                                            Push-Open
     uicontrol('Parent',m_fig,'Style','pushbutton','Units','points',...
                   'position',[ps1-2 9 40 15],...
                   'Tag','S_open1',...
                   'String','Open',...
                   'FontName','Arial','FontSize',8,...  
                   'CallBack','pme(''cb_large'')');
                
               
     %                                            Push-Refresh
     uicontrol('Parent',m_fig,'Style','pushbutton','Units','points',...
                   'position',[ps1+0.5*ps3-20 9 40 15],...
                   'Tag','S_refr1',...
                   'string','Refresh',...
                   'FontName','Arial','FontSize',8,...  
                   'CallBack',...
                   ['if exist(''pmexx2.mat'')==2,','load pmexx2; end;',...
                    'save pmexx1;','pme(''cb_refr'')']);
     %                                            Push-EXIT
     uicontrol('Parent',m_fig,'Style','pushbutton','Units','points',...
                   'position',[ps1+ps3-38 9 40 15],...
                   'Tag','S_exit1',...
                   'string','Close',...
                   'FontName','Arial','FontSize',8,...  
                   'CallBack',[...
                      'if exist(''pmexx1.mat'')==2, delete pmexx1.mat; end;',...
                      'if exist(''pmexx2.mat'')==2, load pmexx2;  delete pmexx2.mat; end;',...
                      'pme(''cb_ex'')']);
     m_f_ch=get(m_fig,'Children');
     for ijk=m_f_ch',
        if (ijk~=menu1) & (ijk~=menu2)
           set(ijk,'Units','normalized');
        end;
     end;
     evalin('caller','save pmexx1');
     pme('cb_type1');
     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-MENU1
   case {'cb_type1','cb_type2','cb_type3'},
     menu1=findobj('Tag','S_pop1');
     all_menu1=get(menu1,'Children');
     set(all_menu1','Checked','off');
     small_w_t1=findobj('Tag','S_type1');
     mAt_tYPe=[]; eval(['mAt_tYPe=',command_str(8),';']);
     set(all_menu1(4-mAt_tYPe),'Checked','on');
     set(small_w_t1,'UserData',mAt_tYPe);
     menu2=findobj('Tag','S_var1');
     small_w_l=findobj('Tag','S_list1');
     set(small_w_l,'Value',1);
     user_struct=get(small_w_l,'UserData');
     user_struct(2:end)=[];
     ch_gl=get(menu2,'UserData');
     
     Pos_STrucT=get(findobj('Tag','Small_PmE'),'UserData');
     sTrLG=Pos_STrucT.strlg; sTrx1=Pos_STrucT.s1;
     sTrx2=Pos_STrucT.s2; sTrx3=Pos_STrucT.s3;
     
     StR0=char(32*ones(1,sTrLG));     
     StR1=StR0; StR1(1:3)=['New']; StR1(sTrx1:sTrx1+5)=['3-by-3'];
     StR1(sTrx2:sTrx2+2)=['pol']; StR1(sTrx3:(sTrx3+length(ch_gl)-1))=ch_gl;   
     user_struct(1).type='pol'; user_struct(1).symb=ch_gl; 
     switch mAt_tYPe,
      case 1, set(small_w_t1,'String','POL matrices.'); 
      case 2, set(small_w_t1,'String','DOUBLE matrices.');
              StR1(sTrx2:sTrx2+5)=['double']; StR1(sTrx3:sTrx3+3)='    ';
              user_struct(1).type='double'; user_struct(1).symb='';
              user_struct(1).matr=[];
      case 3, set(small_w_t1,'String','POL & DOUBLE matrices.');
      otherwise,       
     end
     mAt_list=who('-file','pmexx1.mat');
     [m_l_row,m_l_col]=size(mAt_list);
     if ~isempty(mAt_list),
        for mAt_ind=1:m_l_row
            str_aux=['load pmexx1 ',mAt_list{mAt_ind,1}];
            mTrX=[]; eval(str_aux);
            str_aux=['mTrX=',mAt_list{mAt_ind,1},';']; eval(str_aux);
            tYPeM=class(mTrX); test_aux=0; 
            [rmTrX,cmTrX]=size(mTrX);
            if strcmp(tYPeM,'pol')&(mAt_tYPe==1 | mAt_tYPe==3),
               test_aux=2;
            elseif strcmp(tYPeM,'double')&(mAt_tYPe==2' | mAt_tYPe==3),
               if ndims(mTrX)==2, test_aux=1; end;
            end;
            if test_aux
               StR2=StR0;
               if length(mAt_list{mAt_ind,1}) > sTrx1-2,
                  StR_pOm=mAt_list{mAt_ind,1};
                  StR2(1:sTrx1-2)=[StR_pOm(1:sTrx1-3),'~'];
               else
                  StR2(1:length(mAt_list{mAt_ind,1}))=mAt_list{mAt_ind,1};
               end
               StR_pOm=[int2str(rmTrX),'-by-',int2str(cmTrX)];
               StR2(sTrx1:(sTrx1-1+length(StR_pOm)))=StR_pOm;
               StR2(sTrx2:(sTrx2-1+length(tYPeM)))=tYPeM;
               if test_aux==2,
                  if ~isempty(mTrX.var),
                     StR2(sTrx3:(sTrx3-1+length(mTrX.var)))=mTrX.var;
                  end;
               end
               StR1=[StR1;StR2];
               [r_s1,c_s1]=size(StR1);
               user_struct(r_s1).name=mAt_list{mAt_ind,1};
               user_struct(r_s1).size=StR_pOm;
               user_struct(r_s1).type=tYPeM;
               user_struct(r_s1).symb='';
               if test_aux==2, user_struct(r_s1).symb=mTrX.var; end
               user_struct(r_s1).matr=mTrX;
            end %test_aux
        end %for mat_ind
     end % if ~isempty(mat_list)

     set(small_w_l,'String',StR1);
     set(small_w_l,'UserData',user_struct);
     set(small_w_l,'Enable','On');
     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-VAR
   case 'cb_var',  
     mn_symb=get(gcbo,'Label');
     menu2=get(gcbo,'Parent');
     all_menu2=get(menu2,'Children');
     set(all_menu2','Checked','off');
     set(gcbo,'Checked','on');
     set(menu2,'UserData',mn_symb);
     small_w_l=findobj('Tag','S_list1');
     user_struct=get(small_w_l,'UserData');
     pos_struct=get(findobj('Tag','Small_PmE'),'UserData');
     s3=pos_struct.s3;
     if strcmp(get(small_w_l,'Enable'),'on'),
        small_w_t1=findobj('Tag','S_type1');
        if (get(small_w_t1,'UserData'))~=2,
            StR1=get(small_w_l,'String');
            StR1(1,s3:s3+3)='    ';
            StR1(1,s3:(s3-1+length(mn_symb)))=mn_symb;
            set(small_w_l,'String',StR1);
            user_struct(1).symb=mn_symb;
            set(small_w_l,'UserData',user_struct);
        end;
     end;
     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-POP3
   case 'cb_pop',
     small_w_p3=findobj('Tag','S_pop3');
     small_w_p2=findobj('Tag','S_pop2');
     if get(small_w_p3,'Value')==2,
        set(small_w_p2,'Value',7);
     elseif get(small_w_p2,'Value')==7,
        ch_gl=get(findobj('Tag','S_var1'),'UserData');
        val_var_str=['s','d','p','q','z','z^-1'];         
        val_var=findstr(ch_gl,val_var_str); val_var=val_var(1);
        set(small_w_p2,'Value',val_var);
     end;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-REFRESH
   case 'cb_refr',
     small_w_t1=findobj('Tag','S_type1');
     numxx=int2str(get(small_w_t1,'UserData'));
     eval(['pme(''cb_type',numxx,''');']);
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-LIST
   case 'cb_list',
     fl_click=get(gcf,'SelectionType');
     small_w_l=findobj('Tag','S_list1');
     val_list=get(small_w_l,'Value');
     StR1=get(small_w_l,'String');
     StR1=StR1(val_list,:);
     user_struct=get(small_w_l,'UserData');
    
     switch fl_click, 
      case 'open',
         pme('cb_large');         
      otherwise,
         name_aux=user_struct(val_list).name;
         size_aux=user_struct(val_list).size;
         type_aux=user_struct(val_list).type;
         symb_aux=user_struct(val_list).symb;
         matr_aux=user_struct(val_list).matr;
         val_type=2;
         if strcmp(type_aux,'pol'), val_type=1; end
         if ~isempty(symb_aux)
            val_symb=findstr(symb_aux,['s','d','p','q','z','z^-1']);
            val_symb=val_symb(1);
            set(findobj('Tag','S_pop2'),'Value',val_symb);
         else,
            set(findobj('Tag','S_pop2'),'Value',7);     
         end
         set(findobj('Tag','S_edit1'),'String',name_aux);
         set(findobj('Tag','S_edit1'),'UserData',matr_aux);
         set(findobj('Tag','S_edit2'),'ForegroundColor',[0 0 0]);
         set(findobj('Tag','S_edit2'),'String',size_aux);
         set(findobj('Tag','S_pop3'),'Value',val_type);
      end
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-WARNING1
   case 'cb_warn1', 
     small_w_e2=findobj('Tag','S_edit2');
     small_w_op=findobj('Tag','S_open1');
     set(small_w_op,'Enable','off');
     set(small_w_e2,'Enable','off');
     set(small_w_e2,'String','ERROR !');
     pause(1.2);
     set(small_w_e2,'String','Try again!');
     pause(1.2);
     set(small_w_e2,'String','3-by-3');
     set(small_w_e2,'ForegroundColor',[1 0 0]);
     set(small_w_op,'Enable','on');
     set(small_w_e2,'Enable','on');

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-RESIZE2
   case 'cb_rsf2',
     m_fig=gcbo;
     set(m_fig,'units','points');
     mf_size=get(m_fig,'Position'); mf3=mf_size(3); mf4=mf_size(4);
     set(m_fig,'Units','normalized');
     pos_l1=get(findobj('Tag','S_list1'),'Position');
     if ~isempty(pos_l1),
      pos_struct=get(m_fig,'UserData');
      char1=pos_struct.char1; ls=pos_struct.ls;
      zd=pos_struct.zd;
      mfba=6*ls+42;
      if mf4>mfba,
         ps2=31; ps4=mf4-42-4*ls;
         vt2n1=1.5; vt2n2=9; vt2n3=29;  
         vt2n5=ps2+ps4+1; vt2n6=vt2n5+1;
         vt2n7=ps2+ps4+4+ls; vt2n8=mf4-2*ls;
         
         vt4n1=mf4-3; vt4n2=15; vt4n3=mf4-34-2*ls;
         vt4n5=ls+1; vt4n6=ls; 
         vt4n7=ls; vt4n8=ls;
      else,
         mfcf=mf4/mfba;
         ps2=31*mfcf; ps4=mfcf*(mfba-42-4*ls);
         vt2n1=1.5*mfcf; vt2n2=9*mfcf; vt2n3=29*mfcf;
         vt2n5=ps2+ps4+mfcf; vt2n6=vt2n5+mfcf;
         vt2n7=ps2+ps4+mfcf*(4+ls); vt2n8=mfcf*(mfba-2*ls);
         
         vt4n1=mfcf*(mfba-3); vt4n2=15*mfcf; vt4n3=mfcf*(mfba-34-2*ls);
         vt4n5=mfcf*(ls+1); vt4n6=ls*mfcf; 
         vt4n7=vt4n6; vt4n8=vt4n6;
      end;
      mfbb=36*char1+4*zd+2*ls+47;
      if mf3>mfbb,
         strlg=fix((mf3-40-zd-ls)/char1);
         ps1=20; ps3=strlg*char1+zd+ls+1; 
         s3=strlg-3; s2=strlg-9-fix((2*zd+2+ls)/char1);
         s1=strlg-22-fix((3*zd+4+ls)/char1);
         
         ht1n1=3; ht3n1=mf3-6;
         ht1n21=18; ht1n22=ps1+0.5*ps3-20; ht1n23=ps1+ps3-38; ht3n2=40;
         ht1n3=18; ht3n3=ps3+4;
         ht3n61=4*char1+ls+zd; ht1n61=ps1+ps3-1-ht3n61;
         ht3n62=6*char1+ls+zd; ht1n62=ht1n61-2-ht3n62;
         ht3n51=13*char1+zd; ht1n51=ht1n62-2-ht3n51;
         ht3n52=ps3-ht3n61-ht3n62-ht3n51-7; ht1n52=ps1;
         ht1n8=ps1; ht3n8=ps3;
      else,
         mfco=mf3/mfbb;
         strlg=fix((mfbb-40-zd-ls)/char1);
         ps1=20*mfco; ps3=mf3-40*mfco; 
         s3=strlg-3; s2=strlg-9-fix((2*zd+2+ls)/char1);
         s1=strlg-22-fix((3*zd+4+ls)/char1);
          
         ht1n1=3*mfco; ht3n1=mf3-6*mfco;
         ht1n21=18*mfco; ht1n22=ps1+0.5*ps3-20*mfco;
         ht1n23=ps1+ps3-38*mfco; ht3n2=40*mfco;
         ht1n3=18*mfco; ht3n3=ps3+4*mfco;
         ht3n61=(4*char1+ls+zd)*mfco; ht1n61=ps1+ps3-mfco-ht3n61;
         ht3n62=(6*char1+ls+zd)*mfco; ht1n62=ht1n61-2*mfco-ht3n62;
         ht3n51=(13*char1+zd)*mfco; ht1n51=ht1n62-2*mfco-ht3n51;
         ht3n52=ps3-ht3n61-ht3n62-ht3n51-7*mfco; ht1n52=ps1;
         ht1n8=ps1; ht3n8=ps3;
      end;
      pos_struct.strlg=strlg; pos_struct.s1=s1;
      pos_struct.s2=s2;       pos_struct.s3=s3;
      set(m_fig,'UserData',pos_struct);
      menu1=findobj('Tag','S_pop1');
      menu2=findobj('Tag','S_var1');
      m_f_ch=get(m_fig,'Children');
      for ijk=m_f_ch',
        if (ijk~=menu1) & (ijk~=menu2)
           set(ijk,'Units','points');
        end;
      end;
      set(findobj('Tag','S_frfig1'),'Position',[ht1n1,vt2n1,ht3n1,vt4n1]);    
      set(findobj('Tag','S_frta1'),'Position',[ht1n3,vt2n3,ht3n3,vt4n3]);
      set(findobj('Tag','S_name1'),'Position',[ht1n52,vt2n7,ht3n52,vt4n7]);
      set(findobj('Tag','S_size1'),'Position',[ht1n51,vt2n7,ht3n51,vt4n7]);
      set(findobj('Tag','S_class1'),'Position',[ht1n62,vt2n7,ht3n62,vt4n7]);   
      set(findobj('Tag','S_symb2'),'Position',[ht1n61,vt2n7,ht3n61,vt4n7]);
      set(findobj('Tag','S_pop3'),'Position',[ht1n62,vt2n6,ht3n62,vt4n6]);
      set(findobj('Tag','S_pop2'),'Position',[ht1n61,vt2n6,ht3n61,vt4n6]);
      set(findobj('Tag','S_edit1'),'Position',[ht1n52,vt2n5,ht3n52,vt4n5]);
      set(findobj('Tag','S_edit2'),'Position',[ht1n51,vt2n5,ht3n51,vt4n5]);
      set(findobj('Tag','S_list1'),'Position',[ps1, ps2, ps3, ps4]);
      set(findobj('Tag','S_open1'),'Position',[ht1n21,vt2n2,ht3n2,vt4n2]);
      set(findobj('Tag','S_refr1'),'Position',[ht1n22,vt2n2,ht3n2,vt4n2]);
      set(findobj('Tag','S_exit1'),'Position',[ht1n23,vt2n2,ht3n2,vt4n2]);
      set(findobj('Tag','S_type1'),'Position',[ht1n8,vt2n8,ht3n8,vt4n8]);
      
      for ijk=m_f_ch',
        if (ijk~=menu1) & (ijk~=menu2)
           set(ijk,'Units','normalized');
        end;
      end;
      evalin('caller','save pmexx1');
      pme('cb_refr');
     end; 
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-RESIZE1
   case 'cb_rsf1',
     l_fig=gcbo; flG=0;
     set(l_fig,'Units','points'); 
     pp1=get(l_fig,'Position');
     x_pb=40/pp1(3); y_pb=15/pp1(4);
     set(l_fig,'Units','normalized');
     if (0.101+3*x_pb)<(0.9-2*x_pb), flG=1; end;
     b_obj=findobj(l_fig,'Tag','L_edit1');
     if isempty(b_obj),
        b_obj=b_obj(1); b_str=get(b_obj,'String');
        set(b_obj,'String',b_str);
     end;
     s_figs=get(l_fig,'children');
     for fgr=s_figs'
       if flG & strcmp(get(fgr,'Style'),'pushbutton'),
         switch get(fgr,'String'),
            case 'Browse',
               set(fgr,'Position',[0.9-2*x_pb 0.03 x_pb y_pb]);
            case 'Save as',
               set(fgr,'Position',[0.1+x_pb 0.03 x_pb y_pb]);
            case 'Save', 
               set(fgr,'Position',[0.05 0.03 x_pb y_pb]);
            case 'Close',
               set(fgr,'Position',[0.95-x_pb 0.03 x_pb y_pb]);
          end
       elseif flG & strcmp(get(fgr,'Tag'),'L_ed_saveas'), 
          set(fgr,'Position',[0.101+2*x_pb 0.03 x_pb y_pb]);
       elseif strcmp(get(fgr,'Tag'),'L_edit1'),
          par_cL=get(fgr,'Backgroundcolor'); 
          par_cL=[par_cL(2),par_cL(2),par_cL(3)];
          set(fgr,'Backgroundcolor',par_cL);
          par_eT=get(fgr,'Extent');
          par_pS=get(fgr,'Position');
          if 1.1*par_eT(3)>par_pS(3) | 1.1*par_eT(4)>par_pS(4),
             set(fgr,'Backgroundcolor',(par_cL+[0.15 0 0]));
          end;
          
       end
     end;
  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-LARGE
   case 'cb_large', 
     tag_xx=get(gcbo,'Tag');
     set(findobj('Tag','S_edit2'),'ForegroundColor',[0 0 0]);
     if strcmp(tag_xx,'S_list1'),
       small_w_l=findobj('Tag','S_list1');
       val_list=get(small_w_l,'Value');
       user_struct=get(small_w_l,'UserData');
       Aux_struct=user_struct(val_list);
       Mrx_old=Aux_struct.matr;
     else,
       Aux_struct.name=get(findobj('Tag','S_edit1'),'String');
       Aux_struct.size=get(findobj('Tag','S_edit2'),'String');
       Aux_struct.type='pol';
       if (get(findobj('Tag','S_pop3'),'Value'))==2,
          Aux_struct.type='double';
       end;
       small_w_p2=findobj('Tag','S_pop2');
       vrx_str=get(small_w_p2,'String');
       vrx_val=get(small_w_p2,'Value');
       Aux_struct.symb=deblank(vrx_str(vrx_val,:));
       Mrx_old=get(findobj('Tag','S_edit1'),'UserData');
       Aux_struct.matr=Mrx_old;
     end;
     l_figs=get(0,'children');
     f_lar_mat=[];
     for lf=l_figs',
        if strcmp(get(lf,'Tag'),'Large_mat'),
           f_lar_mat=[f_lar_mat,lf];
        end
     end;
     set(findobj('Tag','S_size1'),'UserData',[]);
     symb_xx=Aux_struct.symb;
     if strcmp(class(Mrx_old),'pol'),
        symb_xx=[symb_xx,Mrx_old.var,get(findobj('Tag','S_var1'),'UserData')];
     else,
        symb_xx=[symb_xx,get(findobj('Tag','S_var1'),'UserData')];   
     end;
     symb_xx=strrep(symb_xx,'z^-1','1'); symb_xx=symb_xx(1);
     Aux_struct.symb=strrep(symb_xx,'1','z^-1');
     Mrx_str=Aux_struct.size;
     Mrx_str=strrep(Mrx_str,'-by-',' ');
     [Mrx_StR1,rem_StR1]=strtok(Mrx_str);
     [Mrx_StR2,rem_StR2]=strtok(rem_StR1);
     Mrx_str=[Mrx_StR1,',',Mrx_StR2];
     Mrx_new=[];
     au_str=['set(findobj(''Tag'',''S_size1''),''UserData'',1); '];
     evalin('caller',...
        ['set(findobj(''Tag'',''S_edit2''),''UserData'',zeros(',Mrx_str,'));'],...
        [au_str,'pme(''cb_warn1'');']);
     if isempty(get(findobj('Tag','S_size1'),'UserData')),  
      Mrx_new=get(findobj('Tag','S_edit2'),'UserData');
      m_r=size(Mrx_new,1);   m_c=size(Mrx_new,2);  
      if m_r>30 | m_c>35,
        m_r2=min(m_r,30); m_c2=min(m_c,35); 
        n_war_small=findobj('Tag','Small_PmE');
        n_war2=figure('numbertitle','Off','Name','W A R N I N G',...
           'Units','normalized','MenuBar','none',...
           'Position',[0.38 0.65 0.4 0.18],'WindowStyle','modal');
        uicontrol(n_war2,'Style','frame','Units','normalized',...
           'Position',[0.01 0.01 0.98 0.98]);
        war_str_n2=char('Matrix is too large and will be cut',...
          ['to first  ',int2str(m_r2),'  rows and  ',int2str(m_c2),'  columns!']);
        uicontrol(n_war2,'Style','text','Units','normalized',...
           'BackgroundColor',[0.5 0.55 0.7],...
           'ForegroundColor',[1 1 1],...
           'FontSize',10,'FontWeight','bold',...
           'Position',[0.1 0.47 0.8 0.32],'String',war_str_n2);
        uicontrol(n_war2,'Style','pushbutton','Units','normalized',...
           'BackgroundColor',[0.5 0.55 0.7],...
           'ForegroundColor',[1 1 1],...
           'UserData',[int2str(m_r2),'-by-',int2str(m_c2)],...
           'String','OK','Position',[0.1 0.1 0.3 0.2],...
           'CallBack',...
           ['set(findobj(''Tag'',''S_edit2''),''String'',get(gcbo,''UserData''));',...
            'close(gcf); pme(''cb_large'')']);
         uicontrol(n_war2,'Style','pushbutton','Units','normalized',...
           'BackgroundColor',[0.5 0.55 0.7],...
           'ForegroundColor',[1 1 1],...
           'String','CANCEL','Position',[0.6 0.1 0.3 0.2],...
           'CallBack',['close(gcf); figure(findobj(''Tag'',''Small_PmE''));']);
         return;
      end;
      Mrx_new=zeros(m_r,m_c);
      Aux_struct.size=[int2str(m_r),'-by-',int2str(m_c)];
      set(findobj('Tag','S_edit2'),'String',Aux_struct.size);
      if length(f_lar_mat) < 4,
       if m_c==0 | m_r==0,
         xx_pos=0.11-0.02*length(f_lar_mat); xx_lg=0.48;
         yy_pos=0.7-0.03*length(f_lar_mat); yy_lg=0.2;
       else,
         if m_c<4, xx_pos=0.11-0.02*length(f_lar_mat); xx_lg=0.48;
         elseif m_c>7, xx_pos=0.11-0.02*length(f_lar_mat); xx_lg=0.84;
         else, xx_pos=0.11-0.02*length(f_lar_mat); xx_lg=0.12*m_c;
         end;
         if m_r<2, yy_pos=0.7-0.03*length(f_lar_mat); yy_lg=0.2;
         elseif m_r>10, yy_pos=0.25-0.03*length(f_lar_mat); yy_lg=0.65;
         else, yy_pos=0.75-0.03*length(f_lar_mat)-0.05*m_r;
               yy_lg=0.15+0.05*m_r;
         end;
       end;
       set(findobj('Tag','Small_PmE'),'Pointer','watch');
       l_fig=figure('numbertitle','Off','units','normalized',...
        'Position',[xx_pos yy_pos xx_lg yy_lg],...
        'ResizeFcn','pme(''cb_rsf1'');',...  
        'name',['Matrix  ',Aux_struct.name],'MenuBar','none',...
        'Tag','Large_mat','UserData',Aux_struct);
       uicontrol(l_fig,'Style','frame','units','normalized',...
        'BackgroundColor',[0.45 0.5 0.65],...
        'ForegroundColor',[0 0 0],...
        'position',[0.01 0.01 0.98 0.98]);
       set(l_fig,'Units','points'); 
       pp1=get(l_fig,'Position');
       x_pb=40/pp1(3); y_pb=15/pp1(4);
       set(l_fig,'Units','normalized');
       %                                             TEXT-StatusWord
       uicontrol(l_fig,'Style','text','units','normalized',...
        'Position',[0.05 0.04+y_pb 0.88 0.149],...
        'BackgroundColor',[0.45 0.5 0.65],...
        'ForegroundColor',[1 1 1],...
        'Tag','Status_wd','String',' ');
       %                                             PUSH-Browse
       uicontrol(l_fig,'Style','pushbutton','units','normalized',...
        'Position',[0.9-2*x_pb 0.03 x_pb y_pb],...
        'FontName','Arial','FontSize',8,...  
        'String','Browse','CallBack',[...
        'figure(findobj(''Tag'',''Small_PmE''));',...  
        'if exist(''pmexx2.mat'')==2, load pmexx2; end;',...
        'save pmexx1; pme(''cb_refr'')']);
       %                                             PUSH-Save as
       uicontrol(l_fig,'Style','pushbutton','units','normalized',...
        'Position',[0.1+x_pb 0.03 x_pb y_pb],'Tag','L_saveas',...
        'String','Save as',...
        'FontName','Arial','FontSize',8,...  
        'CallBack','pme(''cb_saveas''); ');
       %                                             EDIT-Save as
       uicontrol(l_fig,'Style','edit','units','normalized',...
        'Position',[0.101+2*x_pb 0.03 x_pb y_pb],'Tag','L_ed_saveas',...
        'BackgroundColor',[0.9 0.9 1],'ForegroundColor',[0 0 0],...          
        'Enable','off','Visible','off','String',' ',...
        'UserData',[],...
        'CallBack',[...
        'if exist(''pmexx2.mat'')==2, load pmexx2; end; save pmexx1; ',...
        'pme(''cb_saveed''); ',...
        'if exist(''pmexx2.mat'')==2, load pmexx2;  delete pmexx2.mat; end;',...
        'save pmexx1; pme(''cb_refr'')']);
       %                                             PUSH-Save
       uicontrol(l_fig,'Style','pushbutton','units','normalized',...
        'Position',[0.05 0.03 x_pb y_pb],'Tag','L_save',...
        'String','Save',...
        'FontName','Arial','FontSize',8,...  
        'CallBack',['pme(''cb_save''); ',...
        'if exist(''pmexx2.mat'')==2, load pmexx2;  delete pmexx2.mat; end;',...
        'save pmexx1; pme(''cb_refr'')']);
       %                                             PUSH-Close
       uicontrol(l_fig,'Style','pushbutton','units','normalized',...
        'Position',[0.95-x_pb 0.03 x_pb y_pb],...
        'String','Close',...
        'FontName','Arial','FontSize',8,...  
        'CallBack','pme(''cb_close'');');
     
       if strcmp(Aux_struct.type,'pol'),
            l_StR1=[Aux_struct.size,...
                   ' POL matrix in variable ''',Aux_struct.symb,'''.'];
       else,
            l_StR1=[Aux_struct.size,...
                   ' DOUBLE matrix.'];
       end;
       l_str3=int2str(1+length(f_lar_mat));
       L_text_x=uicontrol(l_fig,'Style','text','units','normalized',...
        'BackgroundColor',[0.45 0.5 0.65],...
        'ForegroundColor',[1 1 1],...
        'FontSize',10,...
        'position',[0.1 0.851 0.8 0.138],...
        'Tag',['L_text',l_str3],...
        'UserData',0,...
        'String',l_StR1);
       if m_r==0 | m_c==0,
          uicontrol(l_fig,'Style','text','units','normalized',...
             'BackgroundColor',[0.65 0.7 0.85],...
             'ForegroundColor',[0 0 0],...
             'FontSize',13,...
             'Position',[0.2 0.4 0.6 0.2],...
             'String','Empty Matrix');
           Aux_struct.matr=Mrx_new;        
           set(l_fig,'UserData',Aux_struct);
        
       else,
          I_row=0; J_col=0;
          [m_r_old,m_c_old]=size(Mrx_old);
          Mrx_old_t=class(Mrx_old); Mrx_new_t=Aux_struct.type;
          I_row=1:min(m_r,m_r_old);
          J_col=1:min(m_c,m_c_old);
          Mrx_new=pol(Mrx_new);
          if ~isempty(I_row)&(~isempty(J_col)),
             if strcmp(Mrx_new_t,'double')& strcmp(Mrx_old_t,'pol'),
                Mrx_old=Mrx_old{0};
             end;   
             [wrx1,wrx2]=warning; warning off;
             Mrx_new(I_row,J_col)=Mrx_old(I_row,J_col);
             eval(['warning ',wrx1,';']);
          end;
          pprop(Mrx_new,Aux_struct.symb);
          Aux_struct.matr=Mrx_new;
          set(l_fig,'UserData',Aux_struct);
          C_new=char(Mrx_new);       
          x_width=0.85/(m_c); y_width=(0.61-y_pb)/(m_r);
          c_lg=max(size(Mrx_new));
          for ii=1:m_r
             for jj=1:m_c
                x_dist=0.05+((jj-1)*0.9/(m_c));
                y_dist=0.85-((ii)*(0.66-y_pb)/(m_r));
                if ii<=max(I_row)&jj<=max(J_col),
                   B_color=[0.85 0.85 0.95];
                else,
                   B_color=[0.75 0.75 0.85];
                end;
                if c_lg>1, c_str=C_new{ii,jj};
                else, c_str=C_new; end;
                dat_struct=struct('pos',[ii,jj],'data',Mrx_new(ii,jj),...
                           'sz',0,'hndl',0);
                eDWd=uicontrol('Parent',l_fig,'units','normalized',...
                   'Backgroundcolor',B_color,...
                   'Foregroundcolor',[0 0 0],...
                   'ButtonDownFcn','pme(''cb_down1'');',...   
                   'Enable','inactive',...                 
                   'Position',[x_dist y_dist x_width y_width],...
                   'String',c_str,...
                   'Style','edit',...
                   'UserData',dat_struct,...
                   'Tag','L_edit1');
                sTr_W_h=get(eDWd,'Extent');
                if (1.1*sTr_W_h(3)>x_width)|(1.1*sTr_W_h(4)>y_width),
                   set(eDWd,'Backgroundcolor',(B_color+[0.15 0 0]));
                   dat_struct.sz=1;
                end;
                dat_struct.hndl=eDWd;
                set(eDWd,'UserData',dat_struct);
               end; %for jj=1:m_c
          end; %for ii=1:m_r
       end; %if m_r==0|m_c==0
       set(findobj('Tag','Small_PmE'),'Pointer','arrow');
       
      else,
        l_numb=length(f_lar_mat);
        l_war=figure('numbertitle','Off','Name','W A R N I N G',...
           'Units','normalized','position',[0.25 0.65 0.55 0.18],...
           'Tag','L_warning','MenuBar','none','WindowStyle','modal');
        uicontrol(l_war,'Style','frame','units','normalized',...
           'position',[0.01 0.01 0.98 0.98]);
        war_str=char('To open a new  matrix pad,',...
                     'close one of the existing.');
        uicontrol(l_war,'Style','text','units','normalized',...
           'BackgroundColor',[0.5 0.55 0.7],...
           'ForegroundColor',[1 1 1],...
           'FontSize',10,'FontWeight','bold',...
           'Position',[0.1 0.5 0.8 0.35],...
           'String',war_str,'Tag','L_war_text',...
           'UserData',f_lar_mat);
        for ii=1:l_numb,
           war_str_name=get(f_lar_mat(ii),'Name');
           uicontrol(l_war,'Style','checkbox','Units','normalized',...
              'Position',[0.1+(ii-1)*(0.77/(l_numb)+0.01) 0.27 (0.77/(l_numb)) 0.17],...
              'String',war_str_name,'Value',0,'Tag',['W_',int2str(ii)]);
        end;
        uicontrol(l_war,'Style','pushbutton','Units','normalized',...
           'Position',[0.1 0.07 0.3 0.17],...
           'BackgroundColor',[0.5 0.55 0.7],...
           'ForegroundColor',[1 1 1],...
           'String','O K',...
           'CallBack','save pmexx1; pme(''cb_warn'')');
        uicontrol(l_war,'Style','pushbutton','Units','normalized',...
           'Position',[0.6 0.07 0.3 0.17],...
           'BackgroundColor',[0.5 0.55 0.7],...
           'ForegroundColor',[1 1 1],...
           'String','CANCEL','CallBack','close(gcf);');
      end; %if length(f_lar_mat)<4
     end;
     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-DOWN  
   case 'cb_down1',
     [TeWd,l_fig]=gcbo;
     dat_struct=get(gcbo,'UserData');
     pos_En=dat_struct.pos;
     scr_Uni=get(0,'units'); set(0,'units','centimeters');
     set(TeWd,'units','centimeters'); set(l_fig,'units','centimeters');
     pos_ScR=get(0,'ScreenSize');  set(0,'units',scr_Uni);
     pos_LWd=get(l_fig,'Position');
     pos_eDWd=get(gcbo,'Position');
     sTr_W_h=get(TeWd,'Extent');
     posit_New=pos_eDWd;
     set(TeWd,'units','normalized'); set(l_fig,'units','normalized');
     if (1.1*sTr_W_h(3))>posit_New(3),posit_New(3)=1.1*sTr_W_h(3); end;
     if sTr_W_h(4)>posit_New(4),posit_New(4)=sTr_W_h(4); end;
     posit_New(1)=posit_New(1)+pos_LWd(1);
     posit_New(2)=posit_New(2)+pos_LWd(2);
     posit_New(3)=max(1.5*posit_New(3),4);
     posit_New(4)=max(1.15*posit_New(4),0.85);
     if (posit_New(1)+posit_New(3))>pos_ScR(3),
        posit_New(1)=posit_New(1)+pos_eDWd(3)-posit_New(3);
     end;
     c_str=get(gcbo,'String');
     F_eDWd=figure('units','centimeters',...
                   'Position',posit_New,...
                   'CloseRequestFcn','pme(''cb_clrqf'');',...
                   'Name',['Entry (',int2str(pos_En(1)),',',int2str(pos_En(2)),')'],...
                   'WindowStyle','modal','MenuBar','none',...
                   'UserData',l_fig,...
                   'Tag','L_fig2','NumberTitle','off');
     set(F_eDWd,'units','normalized');
     eDWd=uicontrol(F_eDWd,'Style','edit','units','normalized',...
                   'Backgroundcolor',[1 1 1],...
                   'Foregroundcolor',[0 0 0],...
                   'Position',[0.05 0.05 0.9 0.9],...
                   'String',c_str,...
                   'UserData',dat_struct,...
                   'Tag','L_edit2',...
                   'CallBack','pme(''cb_entr'')');

     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-CLOSEREQUESTF..  
   case 'cb_clrqf',
     F_eDWd=findobj('Tag','L_fig2');
     l_fig=get(F_eDWd,'UserData'); 
     delete(F_eDWd);
     evalin('caller','clear u_dat_xXx_3841 pom_xXxxx_3841;');
     figure(l_fig);
      
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-EXIT  
   case 'cb_ex',
     l_fig=get(0,'children');
     for fgr=l_fig',
        if strcmp(get(fgr,'Tag'),'Large_mat')|strcmp(get(fgr,'Tag'),'Small_PmE'),
           close(fgr);
        end;
     end;
     return;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-WARNING
                                                            %some windows are closed
   case 'cb_warn',
      [h_obj,l_war]=gcbo;
      f_lar_mat=get(findobj(l_war,'Tag','L_war_text'),'UserData');
      l_numb=length(f_lar_mat);
      find_cl=[];
      for ii=1:l_numb,
         if get(findobj(l_war,'Tag',['W_',int2str(ii)]),'Value'),
            find_cl=[find_cl,ii];
         end;
      end;
      if ~isempty(find_cl),
         for ii=find_cl,
            pos_ii=get(f_lar_mat(ii),'Position');
            for jj=1:l_numb,
               pos_jj=get(f_lar_mat(jj),'Position');
               if pos_jj(1)<pos_ii(1),
                pos_new=[pos_jj(1)+0.02, pos_jj(2)+0.03, pos_jj(3), pos_jj(4)];
                set(f_lar_mat(jj),'Position',pos_new);
               end;
            end;
         end
         f_lar_mat=f_lar_mat(find_cl);
         for fg=f_lar_mat,
            close(fg);
         end;
      end;
      close(l_war);
      figure(findobj('Tag','Small_PmE'));  
      pme('cb_refr');
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-CLOSE
   case 'cb_close',
     [h_obj,l_war]=gcbo;                                                       
     l_figs=get(0,'children');
     f_lar_mat=[];
     for lf=l_figs',
        if strcmp(get(lf,'Tag'),'Large_mat'),
           f_lar_mat=[f_lar_mat,lf];
        end
     end;
     find_cl=find(f_lar_mat==l_war);
     pos_ii=get(f_lar_mat(find_cl),'Position');
     for jj=1:(length(f_lar_mat)),
        pos_jj=get(f_lar_mat(jj),'Position');
        if pos_jj(1)<pos_ii(1),
           pos_new=[pos_jj(1)+0.02, pos_jj(2)+0.03, pos_jj(3), pos_jj(4)];
           set(f_lar_mat(jj),'Position',pos_new);
        end;
     end;
     close(l_war); 
     figure(findobj('Tag','Small_PmE'));  
     pme('cb_refr');
                                                            
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-SET NEW
                                                            %         ENTRIES
   case 'cb_entr',
     [ed_wd,F_eDWd]=gcbo;
     a_s=get(ed_wd,'String');
     dat_struct=get(ed_wd,'UserData');
     TeWd=dat_struct.hndl;
     l_fig=get(TeWd,'Parent');
     col_B_g=get(TeWd,'BackgroundColor');
     entr_pos=dat_struct.pos;
     Aux_struct=get(l_fig,'UserData');
     a_s=strrep(a_s,'z^-','@');
     var_xy=['s','d','p','q','z','@'];
     if length(a_s)>1,
        for jj=1:6,
           var_xx=var_xy(jj);
           a_s=strrep(a_s,var_xx,['*',var_xx]);
           b_s=['+';'-';'*';'/';' '];
           for ii=1:5,
              a_s=strrep(a_s,[b_s(ii),'*',var_xx],[b_s(ii),var_xx]);
           end
        end
     end;
     a_s=strrep(a_s,'@','zi^');
     evalin('caller',...
        ['pom_xXxxx_3841=',a_s,';','u_dat_xXx_3841=get(gcbo,''UserData'');',...
           'u_dat_xXx_3841.data=pol(pom_xXxxx_3841);',...
           'set(gcbo,''UserData'',u_dat_xXx_3841);',...
           'clear u_dat_xXx_3841 pom_xXxxx_3841;'],...
        'set(gcbo,''UserData'',[]);');
     str_status=' ';
     dat_struct1=get(ed_wd,'UserData');
     if isempty(dat_struct1), 
        str_status=['Error in new entry definition   (',int2str(entr_pos(1)),...
              ',',int2str(entr_pos(2)),').'];     
        set(ed_wd,'UserData',dat_struct);
     else,
        new_entr=dat_struct1.data;
        new_entr_size=new_entr.size;
        if new_entr_size(1)~=1|new_entr_size(2)~=1,
           str_status=['Error in new entry definition   (',int2str(entr_pos(1)),...
              ',',int2str(entr_pos(2)),').'];     
           set(ed_wd,'UserData',dat_struct);
        end;      
        if strcmp(Aux_struct.type,'double'),
           if new_entr.deg>0,
              str_status=['New entry  (',int2str(entr_pos(1)),',',...
                    int2str(entr_pos(2)),')  is not a constant.'];
              set(ed_wd,'UserData',dat_struct);
           end;
        end;
     end;
     set(findobj(l_fig,'Tag','Status_wd'),'String',str_status);
     Matr_new=Aux_struct.matr;
     vec_col=[1 0 0];
     if strcmp(str_status,' '),
        [wrx1,wrx2]=warning; warning off;
        Matr_new(entr_pos(1),entr_pos(2))=new_entr;
        eval(['warning ',wrx1,';']);
        Aux_struct.matr=Matr_new;
        vec_col=[0 0 0];
        set(l_fig,'UserData',Aux_struct);
     end
     set(ed_wd,'String',char(Matr_new(entr_pos(1),entr_pos(2))));
     set(TeWd,'String',char(Matr_new(entr_pos(1),entr_pos(2))));
     set(ed_wd,'ForegroundColor',vec_col);
     dat_struct.sz=0;
     pOs_W_h=get(TeWd,'Position');
     sTr_W_h=get(TeWd,'Extent');
     set(TeWd,'BackgroundColor',[col_B_g(2),col_B_g(2),col_B_g(3)]);
     if (1.1*sTr_W_h(3)>pOs_W_h(3))|(1.1*sTr_W_h(4)>pOs_W_h(4)),
        set(TeWd,'BackgroundColor',[0.15+col_B_g(2),col_B_g(2),col_B_g(3)]);
        dat_struct.sz=1;
     end;   
     set(ed_wd,'UserData',dat_struct);
     set(TeWd,'UserData',get(ed_wd,'UserData'));
     if strcmp(str_status,' '),
        close(F_eDWd);
        figure(l_fig);
     end;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-SAVE
   case {'cb_save','cb_saveed'},
     [ans,l_fig]=gcbo;
     figure(l_fig);
     New_struct=get(l_fig,'UserData');
     Matr=New_struct.matr;
     if isempty(Matr),
        if(strcmp(New_struct.type,'pol')),  Matr=pol(Matr); end;
     elseif strcmp(New_struct.type,'double'),
        Matr=Matr{:};
     end;
     x=0;
     L_ed_s=findobj(l_fig,'Tag','L_ed_saveas');
     if strcmp(command_str,'cb_save'),
        Name_str=New_struct.name;
     else,
        Name_str=get(findobj(l_fig,'Tag','L_ed_saveas'),'String');
        Name_str=strrep(Name_str,' ','');
        if isempty(Name_str),
           set(L_ed_s,'Enable','off');
           set(L_ed_s,'String','Error!'); pause(1);
           set(L_ed_s,'String','Empty!'); pause(1);
           set(L_ed_s,'String',' ');
           set(L_ed_s,'Enable','on'); return;
        end
        x=pexist(Name_str,'pmexx1.mat');
        ed_struct=struct('matr',Matr,'name',Name_str);
        set(L_ed_s,'UserData',ed_struct);
     end
     if ~x,
        set(findobj(l_fig,'Tag','Status_wd'),'String',' ');
        set(findobj(l_fig,'Tag','Status_wd'),'FontSize',8);
        set(l_fig,'name',['Matrix  ',Name_str]);
        set(L_ed_s,'Enable','off');
        set(L_ed_s,'Visible','off');
        set(findobj(l_fig,'Tag','L_saveas'),'Enable','On');
        psave(Matr,Name_str,'pmexx2.mat');
     else,
        n_war=figure('numbertitle','Off','Name','W A R N I N G',...
           'Units','normalized','MenuBar','none',...
           'Position',[0.38 0.65 0.4 0.18],'WindowStyle','modal');
        uicontrol(n_war,'Style','frame','Units','normalized',...
           'Position',[0.01 0.01 0.98 0.98]);
        war_str_n=char('Matrix name used already exists.',...
                       '        Overwrite?');
        uicontrol(n_war,'Style','text','Units','normalized',...
           'BackgroundColor',[0.5 0.55 0.7],...
           'ForegroundColor',[1 1 1],...
           'FontSize',10,'FontWeight','bold',...
           'Position',[0.1 0.47 0.8 0.32],'String',war_str_n);
        uicontrol(n_war,'Style','pushbutton','Units','normalized',...
           'BackgroundColor',[0.5 0.55 0.7],...
           'ForegroundColor',[1 1 1],...
           'UserData',[n_war,l_fig],...
           'String','YES','Position',[0.1 0.1 0.3 0.2],...
           'CallBack','pme(''cb_over'')');
        str_cb1=['ans=get(gcbo,''UserData''); close(gcf); figure(ans);',...
           'clear ans;',...
           'set(findobj(gcf,''Tag'',''Status_wd''),''String'','' '');',...
           'set(findobj(gcf,''Tag'',''Status_wd''),''FontSize'',8);',...
           'set(findobj(gcf,''Tag'',''L_ed_saveas''),''Enable'',''off'');',...
           'set(findobj(gcf,''Tag'',''L_ed_saveas''),''Visible'',''off'');',...
           'set(findobj(gcf,''Tag'',''L_saveas''),''Enable'',''on'');'];
        uicontrol(n_war,'Style','pushbutton','Units','normalized',...
           'BackgroundColor',[0.5 0.55 0.7],...
           'ForegroundColor',[1 1 1],...
           'UserData',l_fig,...
           'String','CANCEL','Position',[0.6 0.1 0.3 0.2],...
           'CallBack',str_cb1);
     end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-OVERWRITE  
   case 'cb_over',
     figs=get(gcbo,'UserData');
     close(figs(1)); l_fig=figs(2);
     ed_struct=get(findobj(l_fig,'Tag','L_ed_saveas'),'UserData');
     if isstruct(ed_struct),
        psave(ed_struct.matr,ed_struct.name,'pmexx2.mat');
     end;
     figure(l_fig);
     set(findobj(gcf,'Tag','Status_wd'),'String',' ');
     set(findobj(gcf,'Tag','Status_wd'),'FontSize',8);
     set(findobj(gcf,'Tag','L_ed_saveas'),'Enable','off');
     set(findobj(gcf,'Tag','L_ed_saveas'),'Visible','off');
     set(findobj(gcf,'Tag','L_saveas'),'Enable','on');
     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %CALLBACK-SAVE AS     
   case 'cb_saveas',
     [l_saveas,l_fig]=gcbo;
     l_ed=findobj(l_fig,'Tag','L_ed_saveas'); 
     set(l_saveas,'Enable','Off');
     set(l_ed,'String',' ');
     set(l_ed,'Enable','on');
     set(l_ed,'Visible','on');
     
   otherwise,
     error('Invalid input argument.');
    
end; %switch command_str
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% SUBFUNCTION - (pexist) %%%%%%%%%%%%%%%%%%%%
function [exi,A]=pexist(var_str,space_str)
% The function checks if variable with name 'var_str' exists
% in workspace defined by space_str.mat. The output argument is:
%     exi = 1,     if var_str is a variable in space_str.mat
%     exi = 0,     if var_str does not exist.
% The output argument A = var_str.
% NOTE: The input arguments var_str and space_str are strings,
% but output argument A is not a string.
test1=0; exi=0; A=[];
if nargin ~=2
   test1=1;
elseif ~ischar(var_str) | ~ischar(space_str)
   test1=1;
elseif ~strcmp('.mat',space_str(end-3:end))
   test1=1;
end

if test1
    error('Illegal type of input argument.');
end
if exist(space_str)==2,
   clear A exi test1;
   str_aux=['load ',space_str,';',...
         'xX_cHecK_Ex98=exist(''',var_str,''');',...
         'if xX_cHecK_Ex98==1,',...
             'A=',var_str,'; exi=1;',...
         'else, exi=0; A=[];',...
         'end'];
   eval(str_aux);
end;

%end ..pexist

%%%%%%%%%%%%% SUBFUNCTION - (psave) %%%%%%%%%%%%%%%%%%%%
function ans=psave(A,var_name,file_name)
% Function save the input argument A with 'var_name'in MAT-file
% 'file_name.mat'.
% If the MAT-file already exists and contains some variables
% the input variable A is added to them. 
% If in 'file_name.mat' already exists variable with the same
% name 'var_name', the argument A overwrites the old existing
% variable and output argument ans = 1. If A does not overwrite
% any variable ans = 0.
test1=0; 
if nargin ~=3
   test1=1;
elseif ~ischar(file_name) | ~ischar(var_name)
   test1=1;
elseif ~strcmp('.mat',file_name(end-3:end))
   test1=1;
end
if test1
    error('Illegal type of input argument.');
end
B=A;
eval([var_name,'=B;']);
ans=exist(file_name);
if exist(file_name)==2,
   ans=pexist_m(var_name,file_name);
   eval(['save psave_mat ',var_name,' ans;',' clear ',var_name,';']);
   clear test1 B A ans var_name;
   eval(['load ',file_name,';']);
   load psave_mat;
   eval(['clear ans file_name; save ',file_name,';']);
   load psave_mat ans;
   delete  psave_mat.mat;
else,
   ans=0;
   eval(['save ',file_name,' ',var_name,';']);
end;

%end .. psave

%end .. pme
