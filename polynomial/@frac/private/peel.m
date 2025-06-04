function s2 = peel(s1)
%PEEL    Peel off the first and the last two lines in string
%            S2 = PEEL(S1)

%   Author:  J.Jezek  1999
%   Copyright(c) 1999 by Polyx, Ltd.

k = findstr(s1,char(10));
s2 = s1(k+1:end-2);

%end .. @frac/private/peel
