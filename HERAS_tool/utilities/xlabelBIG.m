function xl = xlabelBIG(xlabel_string,varargin)
xl=xlabel(xlabel_string,varargin{:});
xl.Interpreter='latex';
xl.FontSize=14;
end