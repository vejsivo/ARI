function pplot(A,fig_num)
%PPLOT  Two-dimensional plot of a polynomial matrix
%
% The commmand
%    PPLOT(A) 
% creates a new figure with a two-dimensional plot of the polynomial 
% matrix A. The commmand
%    PPLOT(A,N) 
% puts the plot into the existing figure with number N. The commmand
%    PPLOT(A,'LAST') 
% puts it into the last created or mouse-activated figure.
%
% The color of each cell corresponds to the degree of the corresponding 
% entry of A. Empty white places correspond to zero entries.
%
% The color map may be change with the help of the menu 'Color'.

%       Author(s):  S. Pejchova 22-6-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 14-Sep-1998 12:02:34   $
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
       name_str=['Plot of Two-Sided Polynomial Matrix ',inputname(1)];
    else,
       eval('A=pol(A);','error(''Invalid input.'');');
       t_o = 0;
       name_str=['Plot of Polynomial Matrix ',inputname(1)];
    end;
    degA = A.d; As = A.s; [rA,cA]=size(A);
    if isempty(degA)
       disp('pplot: Input matrix is empty ');
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
        if strcmp(name_x,'Simple_Plot of Pol_Mat'),
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
             set(gca,'Tag',[tag_str(1:11),num_plot]);
          end;
       end;
    else,
       if old_plot,
          figure(m_fig); clf reset;
          set(gcf,'name',name_str,...
                  'Tag','Simple_Plot of Pol_Mat',...
                  'units','normalized','MenuBar','none','UserData',['cool']);
       else,
          m_fig=figure('name',name_str,...
                  'units','normalized','MenuBar','none',...
                  'Tag','Simple_Plot of Pol_Mat',...
                  'position',[0.28 0.25 0.7 0.65],'UserData',['cool']);
       end;
       num_plot=int2str(m_fig);
       menu1=uimenu(m_fig,'Label','Color');
          uimenu(menu1,'Label','Cool','CallBack','pplot(''cool'');');
          uimenu(menu1,'Label','Spring','CallBack','pplot(''spring'');');
          uimenu(menu1,'Label','Summer','CallBack','pplot(''summer'');');
          uimenu(menu1,'Label','Autumn','CallBack','pplot(''autumn'');');
          uimenu(menu1,'Label','Winter','CallBack','pplot(''winter'');');
          uimenu(menu1,'Label','Copper','CallBack','pplot(''copper'');');
          uimenu(menu1,'Label','Pink','CallBack','pplot(''pink'');');
          uimenu(menu1,'Label','Hot','CallBack','pplot(''hot'');');
          uimenu(menu1,'Label','Gray','CallBack','pplot(''gray'');');
       axes('position',[.1 .12 .84 .03],'Tag',['S_SmallAxes',num_plot],'Visible','off');   
       axes('position',[.1 .29 .84 .65],'Tag',['S_LargeAxes',num_plot],'Visible','off'); 
    end;
    
     if isinf(degA), degA=0;, end
     Xo=[[cA-1:-1:0; cA-1:-1:0]; [cA:-1:1; cA:-1:1]; [cA-1:-1:0]];
     Yo=[zeros(1,cA);ones(2,cA); zeros(2,cA)];
     X=Xo;, Y=Yo;
     if rA > 1
       for i=1:rA-1
         X=[X, Xo];  Y=[Y, i*ones(5,cA)+Yo];
       end
     end
     DEGA =deg(A,'ent');
     Zau=(fliplr(DEGA))';, Zo=(Zau(:))'; I=isinf(Zo);
     if isempty(I)==0
        Zo(:,I)=[]; X(:,I)=[]; Y(:,I)=[];
     end
     if isempty(Zo)
        Z=[];, Zm=0;
     else
        Zm=max(abs(Zo));
        Cr=Zo+1;
     end
     %%%%%%%%%%%%%% PALETTE  %%%%%%%%%%%
     str_col=get(m_fig,'UserData'); numx1='(Zm+1);'; 
     if strcmp(str_col,'pink')|strcmp(str_col,'hot')|...
           strcmp(str_col,'gray'), numx1='(Zm+2);';
     end;
     eval(['map1=',str_col,numx1]); 
     map1_a=permute(map1(1:(Zm+1),:),[3,1,2]);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     xtls='XTickLabel';ytls='YTickLabel'; 
     if ~isempty(X),
        xo=([0:Zm; 1:Zm+1; 1:Zm+1; 0:Zm]);
        yo=[zeros(2,Zm+1); 0.1*ones(2,Zm+1)];
        xto=0:Zm; xtlo=[int2str(t_o)];
        for i=1:Zm,xtlo=[xtlo,'|',int2str(i+t_o)]; end;
        hsa=findobj('Tag',['S_SmallAxes',num_plot]);
        axes(hsa);
        set(gca,'XLim',[0,Zm+1],'XTick',xto+0.5,xtls,xtlo,...
           'YTick',[],ytls,[],'FontSize',8,'Visible','On');   
        text('String','Degrees of leading coefficients','FontSize',8,...
           'Position',[0.3*(Zm+1) -0.21]);
        patch('XData',xo,'YData',yo,'FaceColor','flat','CData',map1_a,...
             'Tag',['S_ColBar',num_plot]);
     end;
     xt=0:cA; xtl=[' ']; yt=0:rA; ytl=[' '];
     for i=1:cA, xtl=[xtl,'|',int2str(i)]; end;
     for i=1:rA, ytl=[ytl,'|',int2str(i)]; end; 
     hla=findobj('Tag',['S_LargeAxes',num_plot]);
     axes(hla);
     set(gca,'Visible','on',...
        'XLim',[0 cA],'XGrid','off','XTick',xt-0.5,xtls,xtl,...
        'YLim',[0 rA],'YGrid','off','YTick',yt-0.5,ytls,ytl,'YDir','reverse');
     [xx1,yy1]=ndgrid(xt,yt);  [yy2,xx2]=ndgrid(yt,xt);
     line(xx1,yy1,'Color','k'); line(xx2,yy2,'Color','k');
     text('String','COLUMNS','FontSize',8,'Position',[0.4*cA 1.12*rA]);
     text('String','ROWS','FontSize',8,'Position',[-0.07*cA 0.6*rA],...
        'Rotation',90);

     if ~isempty(X),
        Cn=zeros(1,length(Cr),3); Cn(:,1:length(Cr),:)=map1_a(:,Cr,:);
        C=repmat(Cn,[5,1,1]);
        patch('XData',X,'YData',Y,'FaceColor','flat',...
              'CData',C,'Tag',['S_BottomPlot',num_plot],'UserData',Cr);
     end;
     hold off
   %                                                     CALLBACK FOR 'COLOR'
case {'cool','spring','summer','autumn','winter','copper','pink','hot','gray'},
     [objx,figx]=gcbo;
     set(figx,'UserData',[command_str]);
     Name1=int2str(figx);
     clb=findobj('Tag',['S_ColBar',Name1]);
     if ~isempty(clb),
       numx3=get(clb(1),'CData'); numx2=size(numx3,1); 
       numx1=['(',int2str(numx2),');']; 
       if strcmp(command_str,'pink')|strcmp(command_str,'hot')|...
           strcmp(command_str,'gray'), numx1=['(',int2str(numx2+1),');']; 
       end;
       eval(['map1=',command_str,numx1]); 
       map1_a=permute(map1(1:(numx2),:),[3,1,2]);
        if strcmp(get(clb,'Visible'),'on'), set(clb,'CData',map1_a); end;
       btpl=findobj('Tag',['S_BottomPlot',Name1]);
       if ~isempty(btpl),
         Cr=get(btpl(end),'UserData');
         Cn=zeros(1,length(Cr),3); Cn(:,1:length(Cr),:)=map1_a(:,Cr,:);
         C=repmat(Cn,[5,1,1]);
         set(btpl(end),'CData',C);
       end;
     end;
otherwise
     error('Invalid input.');
end; %switch

%end ..pplot
