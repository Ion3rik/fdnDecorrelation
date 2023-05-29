% Function for plotting heatmap for maxCorrelation %
% Jon Fagerström %
% 5.7.2021 %
function h = plotHeatMap(matrix, limits)
    h = heatmap(matrix);
    grid on;
    colormap(gray)
    h.ColorLimits = limits; 
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
end