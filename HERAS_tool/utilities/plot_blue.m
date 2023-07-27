function p = plot_blue(x,y,varargin)

p=plot(x,y,varargin{:});
p.Color='blue';
p.LineWidth=2;

end