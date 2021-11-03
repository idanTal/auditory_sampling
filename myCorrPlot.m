function [r, p] = myCorrPlot(X, Y,corrtype, plotFlag)
% calculate and plot correlation 
% inputs:
% X - variable 1 (to be plotted on the X-axis)
% Y - variable 2 
% corrtype - string. Can be either 'pearson' (default) or 'spearman'
% plotFlag - 1 - plot reults; 0 - no plot

% outputs:
% r - correlation coefficient
% p - p-value

if size(X,1) < size(X,2)
    X = X';
end
if size(X,1) ~= size(Y,1)
    Y = Y';
end
if ~exist('corrtype','var')
    corrtype = 'pearson';
end
[rtemp, ptemp] = corr(X,Y,'type',corrtype);

if size(rtemp,2) > 1
    r = rtemp(1,2);
    p = ptemp(1,2);
else
    r = rtemp;
    p = ptemp;
end
if plotFlag
    figure
    scatter(X, Y,30,'MarkerEdgeColor',[0 .5 .5], 'MarkerFaceColor',[0 .7 .7])
    P = polyfit(X,Y,1);
    yfit = P(1)*X+P(2);
    hold on;
    plot(X,yfit,'k');
    % axis tight
    h = legend('sham');
    pos = get(h,'position');
    legend('off');

    text(pos(1), pos(2), ['r = ', num2str(r), ', p=', num2str(p)],'units','normalized')
end
