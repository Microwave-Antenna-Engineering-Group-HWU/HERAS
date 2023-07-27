function p = imagescNonRectangularGrid(x,y,z)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

tol=1e-12;

x_unique=unique(x);
y_unique=unique(y);

index=[~(((x_unique(2:end)-x_unique(1:end-1)).^2)<tol),true];
x_unique=x_unique(index);

index=[~(((y_unique(2:end)-y_unique(1:end-1)).^2)<tol),true];
y_unique=y_unique(index);

dx_min=min(x_unique(2:end)-x_unique(1:end-1));
dy_min=min(y_unique(2:end)-y_unique(1:end-1));
d=min(dx_min,dy_min)/2;

x_unique=unique([x_unique,x_unique+d,x_unique-d]);
y_unique=unique([y_unique,y_unique+d,y_unique-d]);

index=[~(((x_unique(2:end)-x_unique(1:end-1)).^2)<tol),true];
x_unique=x_unique(index);

index=[~(((y_unique(2:end)-y_unique(1:end-1)).^2)<tol),true];
y_unique=y_unique(index);

[X,Y]=meshgrid(x_unique,y_unique);

Z=nan(size(X));

for n=1:length(x)
    [~,index]=min((X-x(n)).^2+(Y-y(n)).^2,[],'all');
    Z(index)=z(n);
end

p=imagesc(x_unique,y_unique,Z,'AlphaData',~isnan(Z));

end