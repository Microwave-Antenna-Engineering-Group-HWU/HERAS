function p = plot_green(x,y,varargin)

p=plot(x,y,varargin{:});
p.Color=[0.4660 0.6740 0.1880];
p.LineWidth=2;

end