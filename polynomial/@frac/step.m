function varargout = step(F,varargin)
% FRAC/STEP   This function is a Polynomial Toolbox wrapper 
%   for the Control Systems Toolbox function STEP. This 
%   function enables the STEP command to be called with 
%   Polynomial Toolbox objects (fractions) as input 
%   arguments directly, eg. STEP(1/s+1).

if nargout,
    if nargin>1,
        [varargout{1:nargout}] = step(ss(F),varargin{:});
    else
        [varargout{1:nargout}] = step(ss(F));
    end;
else
    if nargin>1,
        step(ss(F),varargin{:});
    else
        step(ss(F));
    end;

end;