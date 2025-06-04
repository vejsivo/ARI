function pplot3(A,fig_num)
%PPLOT3  Three-dimensional plot of a polynomial matrix
%
% The commmand
%    PPLOT3(A) 
% creates a new figure with a three-dimensional plot of the polynomial 
% matrix A. The commmand
%    PPLOT3(A,N) 
% puts the plot into the existing figure with number N. The commmand
%    PPLOT3(A,'last') 
% puts it into the last created or mouse-activated figure.
%
% The height of each bar corresponds to the degree of the corresponding 
% entry of A. The color of each bar corresponds to the absolute value of
% the leading coefficient of the corresponding entry of A. Empty white 
% places correspond to zero entries.
%
% The color map may be changed by the menu 'Color'. The button 'Rotate' 
% turns a mouse-based interactive rotation of the plot view off or on.

%       Author(s):  S. Pejchova 8-6-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 14-Sep-1998 13:28:34   $
%       $Revision: 3.0 $  $Date: 09-Aug-2000  S. Pejchova   $

% Effect on other properties:
% Output is a standard Matlab plot.

narginchk(1,2);
% error(nargchk(1,2,nargin));	%REMOVED IN NEW MATLABS
if ~ischar(A),
   command_str='initialize';
else,
   command_str=A;
end
switch command_str
   %                                                     INITIALIZATION
case 'initialize',
    if isa(A,'tsp'),
       t_o = A.o;  A = A.p;
       if t_o>0, A=shift(A,t_o); t_o=0; end;
       name_str=['Plot of Two-Sided Polynomial Matrix ',inputname(1)];
    else,
       eval('A=pol(A);','error(''Invalid input for polynomial matrix plot.'');');
       t_o = 0;
       name_str=['Plot of Polynomial Matrix ',inputname(1)];
    end;
    degA = A.d; As = A.s; [rA,cA]=size(A);
    if isempty(degA)
       disp('pplot3: Input matrix is empty ');
       return
    end  
    m_fig=[]; old_plot=0;
    if nargin==2,
      s_figs=get(0,'children');
      if ~isempty(s_figs)
          if ischar(fig_num),
            if strcmp(fig_num,'last'), m_fig=s_figs(1);
            else, error('Invalid input string.');
            end;
         elseif isa(fig_num,'double') & any(s_figs==fig_num),
            flm=find(s_figs==fig_num); m_fig=s_figs(flm(1));
         end;
      end;
      if ~isempty(m_fig),
        name_x=get(m_fig,'Tag');
        if strcmp(name_x,'3_D_Plot of Pol_Mat'),
           old_plot=1;
        else, old_plot=2;
        end;
      end;
    end;
    if old_plot==1,
       figure(m_fig); num_plot=int2str(m_fig);
       set(m_fig,'name',name_str);
       ax_s=get(m_fig,'Children');
         for ax=ax_s',
             if strcmp(get(ax,'Type'),'axes'),
             axes(ax);, cla;
             set(gca,'Visible','off');
             set(gca,'NextPlot','replace');
             tag_str=get(gca,'Tag');
             set(gca,'Tag',[tag_str(1:9),num_plot]);
          end;
       end;
          
    else,
       if old_plot,
          figure(m_fig); clf reset;
          set(gcf,'name',name_str,...
             'Tag','3_D_Plot of Pol_Mat',...
             'units','normalized','MenuBar','none','UserData',['cool']);
       else,
          m_fig=figure('name',name_str,...
             'units','normalized','MenuBar','none',...
             'Tag','3_D_Plot of Pol_Mat',...
             'position',[0.28 0.25 0.7 0.65],'UserData',['cool']);
       end;
       num_plot=int2str(m_fig);
       uicontrol(m_fig,'style','frame','units','normalized',...
                   'position',[0.798 0.023 0.144 0.054],...
                   'ForegroundColor',[1 1 1]);
       uicontrol(m_fig,'style','radiobutton','units','normalized',...
                   'position',[0.8 0.025 0.14 0.05],...
                   'FontSize',8,'String','ROTATE',...
                   'ForegroundColor',[0 0 0],...
                   'CallBack','pplot3(''rot'');','Tag',num_plot);
       menu1=uimenu(m_fig,'Label','Color');
          uimenu(menu1,'Label','Cool','CallBack','pplot3(''cool'');');
          uimenu(menu1,'Label','Spring','CallBack','pplot3(''spring'');');
          uimenu(menu1,'Label','Summer','CallBack','pplot3(''summer'');');
          uimenu(menu1,'Label','Autumn','CallBack','pplot3(''autumn'');');
          uimenu(menu1,'Label','Winter','CallBack','pplot3(''winter'');');
          uimenu(menu1,'Label','Copper','CallBack','pplot3(''copper'');');
          uimenu(menu1,'Label','Pink','CallBack','pplot3(''pink'');');
          uimenu(menu1,'Label','Hot','CallBack','pplot3(''hot'');');
          uimenu(menu1,'Label','Gray','CallBack','pplot3(''gray'');');
       axes('position',[.1 .16 .84 .03],'Tag',['SmallAxes',num_plot],'Visible','off');   
       axes('position',[.1 .29 .84 .65],'Tag',['LargeAxes',num_plot],'Visible','off');
       view(-20, 80); 
     end; 
    
     if isinf(degA), degA=0;, end
     xt=0:cA; xtl=[' ']; yt=0:rA; ytl=[' ']; zt=0:degA; ztl=[int2str(t_o)];
     xtls='XTickLabel';ytls='YTickLabel'; ztls='ZTickLabel';
     for i=1:cA, xtl=[xtl,'|',int2str(i)]; end;
     for i=1:rA, ytl=[ytl,'|',int2str(i)]; end;
     for i=1:degA, ztl=[ztl,'|',int2str(i+t_o)]; end;
     if (t_o<0) & (degA<(-t_o)), zt=[zt,(-t_o)]; ztl=[ztl,'|0']; end;
        
     Xo=[[cA-1:-1:0; cA-1:-1:0]+0.01; [cA:-1:1; cA:-1:1]-0.01;...
         [cA-1:-1:0]+0.01];
     Yo=[0.01*ones(1,cA);0.99*ones(2,cA); 0.01*ones(1,cA);...
         0.01*ones(1,cA)];
     X=Xo;, Y=Yo;
     if rA > 1
        for i=1:rA-1
            X=[X, Xo];
            Y=[Y, i*ones(5,cA)+Yo];
        end
     end
     [DEGA,LCA]=deg(A,'ent');
     Cau=abs(LCA);
     Cm=(max(max(Cau)));
     Zau=(fliplr(DEGA))';, Zo=(Zau(:))';
     I=isinf(Zo);
     Z=[Zo ; Zo ; Zo ; Zo; Zo];
     Cau=(fliplr(Cau))';, Co=(Cau(:))';
     if Cm, Cr=Co/Cm; else, Cr=Co; end
     Cr=floor(10*Cr)+1;
     Cr(Cr>10)=10;
     if isempty(I)==0
        Zo(:,I)=[]; Z(:,I)=[];
        X(:,I)=[]; Y(:,I)=[];  Cr(I)=[];
     end
     if isempty(X)==0
        X1=X(1,:);, X3=X(3,:); Y1=Y(1,:);,  Y2=Y(2,:);
        XLL=[X1; X1; X1; X1; X1]; YLL=[Y1; Y2; Y2; Y1; Y1];
        XFF=[X1; X3; X3; X1; X1]; YFF=[Y2; Y2; Y2; Y2; Y2];
        XRR=[X3; X3; X3; X3; X3]; YRR=[Y2; Y1; Y1; Y2; Y2];
        XBB=[X3; X1; X1; X3; X3]; YBB=[Y1; Y1; Y1; Y1; Y1];
        ZLL=[zeros(2,length(Zo)); Zo; Zo; zeros(1,length(Zo))];
        ZDD=zeros(5,length(Zo));
        XG=[XBB;XRR;XFF;XLL;X]; YG=[YBB;YRR;YFF;YLL;Y];
        ZG=[ZLL;ZLL;ZLL;ZLL;Z]; CG=Cr;
        zind=find(Zo==0);
        if ~isempty(zind),
           XG(:,zind)=[]; YG(:,zind)=[]; ZG(:,zind)=[]; CG(:,zind)=[];
        end;        
        [rx,cx]=size(XG);
        if ~isempty(CG),
           XG=reshape(XG,[5,rx*cx/5]); YG=reshape(YG,[5,rx*cx/5]);
           ZG=reshape(ZG,[5,rx*cx/5]);        
        end;
     end
     %%%%%%%%%%%%%% PALETTE  %%%%%%%%%%%
     str_col=get(m_fig,'UserData'); numx1='(10);'; 
     if strcmp(str_col,'pink')|strcmp(str_col,'hot')|...
           strcmp(str_col,'gray'), numx1='(12);';
     end;
     eval(['map1=',str_col,numx1]); 
     map1_a=permute(map1(1:10,:),[3,1,2]);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if Cm,
        xo=([0:9; 1:10; 1:10; 0:9])*0.1*Cm;
        yo=[zeros(2,10); 0.1*ones(2,10)];
        xto=[0,xo(2,2:2:10)]; xtlo=['0'];
        for i=2:length(xto),
           xtlo=[xtlo,'|',num2str(xto(i))];
        end;
        hsa=findobj('Tag',['SmallAxes',num_plot]);
        axes(hsa);
        set(gca,'XLim',[0,xto(length(xto))],'XTick',xto,xtls,xtlo,...
           'YTick',[],ytls,[],'FontSize',8,'Visible','On'); 
        text('String','Absolute value of leading coefficients','FontSize',8,...
             'Position',[0.3*Cm -0.18]);
        patch('XData',xo,'YData',yo,'FaceColor','flat','CData',map1_a,'Tag',['ColBar',num_plot]);
     end;
     
     hla=findobj('Tag',['LargeAxes',num_plot]);
     axes(hla);
     set(gca,'Visible','on',...
        'XLim',[0 cA],'XGrid','off','XTick',xt-0.5,xtls,xtl,...
        'YLim',[0 rA],'YGrid','off','YTick',yt-0.5,ytls,ytl,'YDir','reverse',...
        'ZLim',[0 max([1,degA,-t_o])],'ZGrid','off','ZTick',zt,ztls,ztl);
     
     [xx1,yy1]=ndgrid(xt,yt);  [yy2,xx2]=ndgrid(yt,xt);
     line(xx1,yy1,-t_o*ones(size(xx1)),'Color','k','LineStyle',':'); 
     line(xx2,yy2,-t_o*ones(size(xx2)),'Color','k','LineStyle',':');
     line(xx1,yy1,zeros(size(xx1)),'Color','k'); 
     line(xx2,yy2,zeros(size(xx2)),'Color','k');
     
     text('String','COLUMNS','FontSize',8,'Position',[0.4*cA 1.2*rA]);
     text('String','ROWS','FontSize',8,'Position',[-0.2*cA 0.6*rA]);
     text('String','DEGREES','FontSize',8,'Position',[-0.08*cA 0],...
        'Rotation',90);
     if isempty(X)==0,
        Cn=zeros(1,length(Cr),3); Cn(:,1:length(Cr),:)=map1_a(:,Cr,:);
        C=repmat(Cn,[5,1,1]);
        patch('XData',X,'YData',Y,'ZData',ZDD,'FaceColor','flat',...
              'CData',C,'Tag',['BottomPlot',num_plot],'UserData',Cr);
        if ~isempty(CG),
           set(gca,'NextPlot','add');   
           CrG=[CG;CG;CG;CG;CG]; CrG=CrG(:)'; CnG=zeros(1,length(CrG),3);      
           CnG(:,1:length(CrG),:)=map1_a(:,CrG,:); Cg=repmat(CnG,[5,1,1]);
           patch('XData',XG,'YData',YG,'ZData',ZG,'FaceColor','flat',...
                 'CData',Cg,'Tag',['LegoPlot',num_plot],'UserData',CrG);
        end
     end;
     hold off
   %                                                     CALLBACK FOR 'ROTATE'
   case 'rot',
     name1=get(gcbo,'Tag');
     val1=get(gcbo,'Value'); 
     hla=findobj('Tag',['LargeAxes',name1]); 
     axes(hla);
     if val1,
        rotate3d on;
     else,
        rotate3d off;
     end;
   %                                                     CALLBACK FOR 'COLOR'
   case {'cool','spring','summer','autumn','winter','copper','pink','hot','gray'},
     numx1='(10);'; 
     if strcmp(command_str,'pink')|strcmp(command_str,'hot')|...
           strcmp(command_str,'gray'), numx1='(12);';
     end;
     eval(['map1=',command_str,numx1]); 
     map1_a=permute(map1(1:10,:),[3,1,2]);
     [objx,figx]=gcbo;
     set(figx,'UserData',[command_str]);
     Name1=int2str(figx);
     clb=findobj('Tag',['ColBar',Name1]);
     if strcmp(get(clb,'Visible'),'on'), set(clb,'CData',map1_a); end;
     btpl=findobj('Tag',['BottomPlot',Name1]);
     if ~isempty(btpl),
        Cr=get(btpl(end),'UserData');
        Cn=zeros(1,length(Cr),3); Cn(:,1:length(Cr),:)=map1_a(:,Cr,:);
        C=repmat(Cn,[5,1,1]);
        set(btpl(end),'CData',C);
     end;
     lgpl=findobj('Tag',['LegoPlot',Name1]);
     if ~isempty(lgpl),
        CrG=get(lgpl(end),'UserData');
        CnG=zeros(1,length(CrG),3);
        CnG(:,1:length(CrG),:)=map1_a(:,CrG,:); Cg=repmat(CnG,[5,1,1]);
        set(lgpl(end),'CData',Cg);
     end;
   otherwise
     error('Invalid input for polynomial matrix plot.');
  end; %switch
  
  %end .. pplot3
  