function [varargout] = find_closest(X,Y,x0,y0)
% This function finds the index of the point in (X,Y) closest to (x0,y0).
% The distance corresponds to the Euclidean distance.
% INPUTS:
% X, Y:     1D or 2D arrays with the same size. The couple (X(i),Y(i))
%           identifies one point in the x-y plane
% x0, y0:   (x0,y0) target point
% OUTPUTS:
% index_array:  index_array corresponds to the index of the element in 
%               (X,Y) closest to (x0,y0).

[~,linear_index]=min((X-x0).^2+(Y-y0).^2,[],'all');
if numel(X)==length(X)
    varargout{1}=linear_index;
else
    [varargout{1},varargout{2}]=ind2sub(size(X),linear_index);
end

end