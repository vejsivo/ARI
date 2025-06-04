function varargout = bode(F,varargin)
% FRAC/BODE   This function is a Polynomial Toolbox wrapper 
%   for the Control Systems Toolbox function BODE. This 
%   function enables the BODE command to be called with 
%   Polynomial Toolbox objects (fractions) as input 
%   arguments directly, eg. bode(1/s+1).

if nargout,
    if nargin>1,
        [varargout{1:nargout}] = bode(ss(F),varargin{:});
    else
        [varargout{1:nargout}] = bode(ss(F));
    end;
else
    if nargin>1,
        bode(ss(F),varargin{:});
    else
        bode(ss(F));
    end;

end;