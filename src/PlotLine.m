function XYdata = PlotLine(point, slope, Xrange)
%point = 1x2 float (x,y)
%slope = 1x1 float
%Xrange = Nx1 float vector
%use Y = mX + b
arguments
   point (1,2) {mustBeFloat}
   slope (1,1) {mustBeFloat}
   Xrange (1,:) {mustBeFloat}
end
nXvals = length(Xrange);
XYmat = zeros(nXvals,2);
b = point(2)-(slope*point(1));
for i=1:nXvals
   XYmat(i,1) = Xrange(i);
   XYmat(i,2) = slope*Xrange(i) + b;
end
XYdata.XYmat = XYmat;
XYdata.Xintercept = b;
if slope ~= 0
   yint = -b/slope;
else
   yint = NaN;
end
XYdata.Yintercept = yint;
end

