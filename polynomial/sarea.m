function S = sarea(varargin)
%SAREA  Stability test for a family of polynomials with parametric uncertainties.
%
% The commmand
%    S=SAREA(Q1,...,QM,EXPRESSION-STRING,P0,P1,...,PN [,TOL])
% returns M-dimensional array S with 1-s (0-s) where the particular 
% polynomial in the family appears to be stable (unstable). 
% The uncertainty structure represented symbolically by the string 
%
% EXPRESSION-STRING = 'P0(s) + Expr1*P1(s) + Expr2*P2(s) + ... + ExprN*PN(s)'
%            where Expr1, Expr2, ...,ExprN are any expressions composed 
%            from the variable names used for parameter vectors 
%            Q1, ..., QM and those are consistent with Matlab syntax, such as
%            EXPRESSION-STRING='P0(s) + Q1*P1(s) + Q2*P2(s) + ...' or
%            EXPRESSION-STRING='P0(s) + (Q1*Q2-Q3)*P1(s)' etc.
%
% Q1, Q2, ..., QM are vectors of (real or complex) parameters, existing
%                 in the current workspace.
% P0, P1, ..., PN are polynomials existing in the current workspace.
% Please, note that the first M input arguments Q1, ..., QM must be 
% the parameter vector names rather than the vector values and also the 
% next N+1 input arguments P0,P1,...,PN (following EXPRESSION-STRING) 
% must be the polynomial matrix names rather than the polynomial values.  
%
% The array S is length(Q1)-by-length(Q2)-by-...-by-length(QM).
%
% A tolerance TOL may be specified as an additional input argument.
%
% See also SAREPLOT.

%       Author(s):  S. Pejchova 12-11-99
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 04-Feb-2000 11:13:34   $
%       $Revision: 2.5 $  $Date: 26-Jul-2000 10:43:34   $
%       $Revision: 3.0 $  $Date: 02-Nov-2000  S. Pejchova   $
%                         $Date: 04-Dec-2002  S. Pejchova   $

global PGLOBAL;
eval('PGLOBAL.VERBOSE;', 'painit;');
n_I=nargin; eX_Str=''; n_Q=0; n_P=0; P_iY=[];  N_iY_1=0; p_V='no';

p_V=PGLOBAL.VERBOSE;
TOl = PGLOBAL.ZEROING;
last_one=varargin{n_I};
if isa(last_one,'double') & all(size(last_one)==1),  
   TOl=last_one; n_I=n_I-1; 
end;
if n_I<3, error('Not enough input arguments.');  end;
C_iY_n={'.\';'.\(';'./';'./(';'.*';'.*('};

fl1=find(cellfun('isclass',varargin,'char'));
if length(fl1)~=1, error('Expression string is missing or too many strings.'); end;

eX_Str=varargin{fl1};
eX_Str=strrep(eX_Str,' ','');, eX_Str=strrep(eX_Str,'*','.*');
eX_Str=strrep(eX_Str,'/','./');, eX_Str=strrep(eX_Str,'\','.\');
eX_Str=strrep(eX_Str,'^','.^');, eX_Str=strrep(eX_Str,'..','.');
n_Q=fl1-1; n_P=n_I-fl1;
varargin(fl1)=[]; n_I=n_I-1;
if n_Q<1 |~n_P,
   error('Not enough input arguments.');
end;
if n_Q>1,
   [varargin{1:n_Q}]=ndgrid(varargin{1:n_Q});
end;
sz_mx=size(varargin{1});
N_iY_1=prod(sz_mx);
for i_I_i=1:n_Q,
   N_iY_n=inputname(i_I_i); N_iY_v=varargin{i_I_i}; N_iY_v=N_iY_v(:);
   if isempty(N_iY_n), error('Parameter vector must be a named variable.'); end;
   eval([N_iY_n,'=N_iY_v;'],'error(peel(lasterr));');
end;
N_iY_2=size(varargin{n_Q+1});
for i_I_i=n_Q+1:n_I,
   N_iY_n=inputname(i_I_i+1);
   if isempty(N_iY_n), error('Polynomial input argument must be a named variable.'); end;
   N_iY_v=varargin{i_I_i};
   if any(N_iY_2-size(N_iY_v)), error('Polynomial matrices not of the same dimensions.'); end
   N_iY_v=N_iY_v(:).';
   eval([N_iY_n,'=N_iY_v;'],'error(peel(lasterr));');
   for j_J_j=1:6,
      c_iX_n=C_iY_n{j_J_j};
      eX_Str=strrep(eX_Str,[c_iX_n,N_iY_n],[c_iX_n(2:end),N_iY_n]);
   end;
end;

eX_s1=cumsum((eX_Str=='(')-(eX_Str==')'));
%eX_s2 -> is 1 for + or - not inside parentheses.
eX_s2=((eX_Str=='-') | (eX_Str=='+')) & (eX_s1==0);
f_E_2=find(eX_s2(2:end));

f_E_2=[0,f_E_2];  V_auX=[];  V_auX1=[];

for kk=fliplr(f_E_2),
   str_XA=['V_auX1=',eX_Str(kk+1:end),';'];
   eval(str_XA,'error(peel(lasterr));');
   eX_Str(kk+1:end)=[];
   if size(V_auX1,1)~=N_iY_1, V_auX1=repmat(V_auX1,[N_iY_1,1]); end;
   if size(V_auX1,2)==1, V_auX1=repmat(V_auX1,[1,prod(N_iY_2)]); end;
   if ~isempty(V_auX), V_auX=V_auX+V_auX1; 
   else,  V_auX=V_auX1; 
   end;
end;   
S=zeros(sz_mx); n_p=N_iY_1;
step_ver=fix(n_p/10);
tm0=clock;
for jj=1:n_p,
   S_aux=V_auX(jj,:);
   S_aux=reshape(S_aux,N_iY_2);
%   eval(['S(jj)=isstable(S_aux,TOl);'],['S(jj)=NaN;']);
   try,
       S(jj)=isstable(S_aux,TOl);,
   catch,
       S(jj)=NaN;
   end

   if strcmp(p_V,'yes'),
      j_step=jj/step_ver;
      if jj==1,
         disp(sprintf(' Stability to be tested in %d  points.',n_p));
      elseif fix(j_step)==(j_step),
         if (jj/n_p)< 0.95,
            tm1=etime(clock,tm0)*(n_p-jj)/jj;
            tm_str=[int2str(fix(tm1/60)),' min. ',num2str(mod(tm1,60),3),' sec.'];
            disp([sprintf(' %d %% (%d points) checked. Remaining time ',...
                  round(100*jj/n_p),jj),tm_str]);
         else, 
            disp(' All points checked.');
         end;
            
      end;
   end;
end;

%end .. sarea
