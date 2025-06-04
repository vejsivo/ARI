function n = norm(A,varargin)
% @TSP/NORM temporary macro that computes the norm of a TSP A via
% macro @POL/NORM: N = norm(pol(shift(A,deg(A))))
% By Didier Henrion, September 1, 2000.
%
% Corrected by J.Jezek 26-Jul-2001,  deg replaced by tdeg
%                      17-Oct-2001

if ~isa(A,'tsp'),
    error('Some argument but not 1st is invalidly tsp.');
end;
eval('n = norm(pol(shift(A,tdeg(A))),varargin{:});',...
    'error(peel(lasterr));');

%end .. @tsp/norm

  
