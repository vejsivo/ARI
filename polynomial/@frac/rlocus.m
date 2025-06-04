function varargout = rlocus(F,varargin)
% FRAC/RLOCUS   This function is a Polynomial Toolbox wrapper 
%   for the Control Systems Toolbox function RLOCUS. This 
%   function enables the RLOCUS command to be called with 
%   Polynomial Toolbox objects (fractions) as input 
%   arguments directly, eg. RLOCUS(1/s+1).

if nargout,
    if nargin>1,
        [varargout{1:nargout}] = rlocus(ss(F),varargin{:});
    else
        [varargout{1:nargout}] = rlocus(ss(F));
    end;
else
    if nargin>1,
        rlocus(ss(F),varargin{:});
    else
        rlocus(ss(F));
    end;
end;