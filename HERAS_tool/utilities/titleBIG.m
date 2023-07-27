function t = titleBIG(titleString,varargin)

t=title(titleString,varargin{:});
t.FontSize=20;
t.Interpreter='latex';

end