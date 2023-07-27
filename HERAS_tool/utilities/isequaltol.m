function equality_flag = isequaltol(in1,in2,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if nargin==3
    tolerance=varargin{1};
else
    tolerance=1e-6;
end

in1=in1(:);
in2=in2(:);

greater_than_tolerance_array=abs((in1-in2)./in1)>tolerance;
equality_flag=sum(greater_than_tolerance_array)==0;

end