function p = plot_red(x,y,varargin)

p=plot(x,y,varargin{:});
p.Color='red';
p.LineWidth=2;

end