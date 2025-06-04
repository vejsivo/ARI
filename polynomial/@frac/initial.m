function varargout = initial(F,varargin)
% FRAC/INITIAL   This function is a Polynomial Toolbox wrapper 
%   for the Control Systems Toolbox function INITIAL. This 
%   function enables the INITIAL command to be called with 
%   Polynomial Toolbox objects (fractions) as input 
%   arguments directly, eg. INITIAL(1/s+1).

if nargout,
    if nargin>1,
        [varargout{1:nargout}] = initial(ss(F),varargin{:});
    else
        [varargout{1:nargout}] = initial(ss(F));
    end;
else
    if nargin>1,
        initial(ss(F),varargin{:});
    else
        initial(ss(F));
    end;

end;