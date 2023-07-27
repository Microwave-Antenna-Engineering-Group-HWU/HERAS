function zl = zlabelBIG(zlabel_string,varargin)
zl=zlabel(zlabel_string,varargin{:});
zl.Interpreter='latex';
zl.FontSize=14;
end