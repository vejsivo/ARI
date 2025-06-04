function H = chsamp(G,varargin)
%CHSAMP     Change sampling of constant
%
% For constant G (i.e. the standard Matlab matrix),
% the command  H = CHSAMP(G,T,TN,TAU,TAUN)
% in all nontrivial cases results in error.
%
% This macro exists only for completeness.
% See also FRAC/CHSAMP.

%      Author: J.Jezek, 05-Dec-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 05-Feb-2001 $
%                    $ Date 25-Jul-2002 $

if nargin<1,
   error('Not enough input arguments.');
end;
if ~isa(G,'double') | ndims(G)~=2,
   error('Invalid 1st argument.');
end;

lv = length(varargin); 
taun = 0;
start = 1;
if lv>=1,
   arg = varargin{1};
   if isa(arg,'double') & ndims(arg)==2,
      if ~isempty(arg),
         if length(arg)==1 & isreal(arg) & arg>=0,
         else
            error('Invalid old sampling period.');
         end;
      end;
      start = 2;
      if lv>=2,
         arg = varargin{2};
         if isa(arg,'double') & ndims(arg)==2,
            if ~isempty(arg),
               if length(arg)==1 & isreal(arg) & arg>=0,
               else
                  error('Invalid new sampling period.');
               end;
            end;
         end;
         start = 3;
         if lv>=3,
            arg = varargin{3};
            if isa(arg,'double') & ndims(arg)==2,
               if ~isempty(arg),
                  if length(arg)==1 & isreal(arg) & arg>=0,
                     tau = arg;
                  else
                     error('Invalid old sampling phase.');
                  end;
               end;
               start = 4;
               if lv>=4,
                  arg = varargin{4};
                  if isa(arg,'double') & ndims(arg)==2,
                     if ~isempty(arg),
                        if any(size(arg)==1) & ...
                              isreal(arg) & all(arg>=0),
                           taun = arg;
                        else
                           error('Invalid new sampling phase.');
                        end;
                     end;
                     start = 5;
                  end;
               end;
            end;
         end;
      end;
   end;
end;

ltaun = length(taun); staun = size(taun);
HH = G;
if ltaun==1,
   H = HH;
else
   H = cell(staun);
   for k = 1:ltaun,
      H{k} = HH;
   end;
end;
if isempty(G) | all(all(G==0)),
   return;
end;

error('Invalid 1st argument; has not required zero point.');

%end .. chsamp
