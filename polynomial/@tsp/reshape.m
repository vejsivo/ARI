function Q = reshape(T,varargin)
%RESHAPE        Change size of two-sided polynomial
%
% For two-sided polynomial matrix T, the command  Q = RESHAPE(P,M,N)
% returns a M-by-N two-sided polynomial matrix whose elements are
% taken columnwise form T. An error results if T has not M*N elements.
%
% Q = RESHAPE(P,[M N])  is the same thing.
%
%      Author:  J. Jezek 04-Apr-2000
%      Copyright(c) by Polyx, Ltd.
%      $ Revision $  $ Date 24-May-2000 $
%                    $ Date 31-Oct-2000 $

eval('T = tsp(T);','error(peel(lasterr));');
ni = length(varargin);
PP = 0;
eval('PP = reshape(T.p,varargin{1:ni});','error(peel(lasterr));');
Q = tsp(PP); Q.h = T.h;
Q = shift(Q,T.o);

%end .. @tsp/reshape
