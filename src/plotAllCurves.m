function plotAllCurves(XYY,titlestring)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
t = tiledlayout(8,12,'TileSpacing','none');
title(t,titlestring);
xlabel(t,"Time (h)");
ylabel(t,"Optical Density (600 nm)");
Ydata = XYY(:,2:end);
ymax_val = max(Ydata,[],"all");
for i=2:97                                   %for each data col
   nexttile                                  %add new tile to graph
   plot(XYY(:,1),XYY(:,i),'.k');             %plot in tile
   ylim([0,ymax_val]);
   set(gca,'xtick',[],'ytick',[]);           %plot a tile
end   
end