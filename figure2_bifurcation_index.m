%% Bifurcation index and topographic plots
clear; close all; clc
%% initialize
subnum = [2 3 4 5 6 7 8 9 10 11 13 14 15 17 18 19 20 21];
subName = {'pilot_IT01_first_version', 'ML02', 'ML03', 'MY04', 'ER05', 'CC06', 'LR07', 'IT08', 'SD09', 'AD10', 'SF11', 'NN12', 'CM13', 'EL14','DA15', 'XX16', 'BJ17', 'JS18', 'AS19', 'SY20', 'AG21'};
% add and set paths 
addpath(genpath('D://Matlab//fieldtrip-20161212'))
addpath(genpath('D://auditory_sampling//source'))
dataRoot = sprintf('D://auditory_sampling//data');
cd(dataRoot)
% END intialization

%% figure 2: target locked bifurcation index statistics
% load timeoi and freqoi
load(sprintf('%s//timeoi', dataRoot))
load(sprintf('%s//freqoi', dataRoot))
% define time for the statistical test
toi = [-.5 .5];
toiI = [it_nearest(timeoi,toi(1)), it_nearest(timeoi,toi(2))];
timeoiNew = timeoi(toiI(1):toiI(2));
% load the bifurcation index data for all subjects
load(sprintf('%s//bifurcation_index_all_subs', dataRoot))
% define 3d matrix with all subjects averaged across channels
BIgaAllChn = squeeze(mean(BIall,2));
% what is the statistical test to be used (ttest or signrank)
stat = 'signrank';
[pValBI, statValBI] = it_matrixStat(BIgaAllChn, stat);
% correct using FDR
[~, ~, ~, adjP_BI]=fdr_bh(pValBI);

% interpolate map for visualization
[n, m] = size(adjP_BI);
[x,y] = meshgrid(timeoiNew,freqoi);% low-res grid
[x2,y2] = meshgrid(timeoiNew(1):1/500/5:timeoiNew(end),freqoi(1):.1:freqoi(end));  % high-res grid
adjP_BIinterp = interp2(x,y,adjP_BI, x2,y2, 'cubic'); % interpolate up
adjStat_BIinterp = interp2(x,y,statValBI, x2,y2, 'cubic'); % interpolate up

figure
obj1 = imagesc(x2(1,:), y2(:,1), adjStat_BIinterp);
cLimit = get(gca,'clim');
set(gca, 'ydir', 'normal','clim',cLimit*.85)
hold on 
alphaLevel = double(adjP_BIinterp<.05);
ylabel('Frequency, Hz')
xlabel('Time, sec.')
colormap parula
axis tight
colorbar('location','east')
set(gca,'xlim',[-.5 .5])
title('BIfurcation Index')
hold on 
% detect clusters and calculate correlation for each cluster
% use labeling to compute connected clusters
clustLabelsAlpha = bwlabel(alphaLevel);
% select the 4 largest clusters before the onset of the target (i.e. before
% time 0)
t0 = it_nearest(x2(1,:),0);
clustLabelsAlphaPre = clustLabelsAlpha(:,1:t0);
clusterSize = zeros(1,max(max(clustLabelsAlphaPre)));
for coi = 1 : max(max(clustLabelsAlphaPre))
    clusterSize(coi) = sum(sum(clustLabelsAlphaPre == coi));
end
[~,ind] = sort(clusterSize,'descend');
for coi = ind(1:4)
    clustLabelsAlphaTemp = clustLabelsAlpha;
    clustLabelsAlphaTemp(clustLabelsAlphaTemp ~= coi) = 0;
    contour(x2(1,:), y2(:,1), clustLabelsAlphaTemp,'k')
end

%% find the peaks in time and frequency for topographic plots
% get the mean BI for pre stimulus significant pixels
toiI = it_nearest(timeoiNew,0);
BI_GApre = BIgaAllChn(:,:,1:toiI);
maskBIpre = adjP_BI < .05;
maskBIpre = double(maskBIpre(:,1:toiI));
maskBIpre(maskBIpre == 0) = nan;
BIpre = zeros(1,length(subnum));
for subI = 1 : length(subnum)
    BIpre(subI) = squeeze(mean(mean(squeeze(BI_GApre(subI,:,:)).*maskBIpre)));
end

% detect clusters and calculate correlation for each cluster
% use labeling to compute connected clusters
maskBIpre = adjP_BI < .05;
maskBIpre = double(maskBIpre(:,1:toiI));
clustLabels = bwlabel(maskBIpre);
% define prestimulus clusters of interest
clusterSizePre = zeros(1,max(max(clustLabels)));
for coi = 1 : max(max(clustLabels))
    clusterSizePre(coi) = sum(sum(clustLabels == coi));
end
[~,ind] = sort(clusterSizePre,'descend');
coi = ind(1:4);
save(sprintf('%s//coiMax', dataRoot), 'coi')
statValClust = cell(1,length(coi));
peakFreqClust = zeros(2, length(coi));
peakTClust = zeros(2, length(coi));
for coiI = 1 : length(coi)
    % find indices of cluster
    [row,col] = find(clustLabels == coi(coiI));
    statValClust{coiI} = zeros(size(statValBI));
    for ii = 1 : length(row)
        statValClust{coiI}(row(ii), col(ii)) = statValBI(row(ii), col(ii));
    end
    % find the time and frequency limits of the cluster
    rowIndClust = [min(row) max(row)];
    colIndClust = [min(col) max(col)];
    % get the frequency and time for the cluster
    peakFreqClust(:,coiI) = freqoi(rowIndClust);
    peakTClust(:,coiI) = timeoiNew(colIndClust);
end

%% plot topography at the peak of each cluster
% load the electrodes layout
load(sprintf('%s//layoutGTecNew', dataRoot))
% load the average bifurcation index for all channels (fieldtrip structure)
load(sprintf('%s//dataBifInd', dataRoot))

% set color limits for the figures
bifIndexZlim = squeeze(median(BIgaAllChn,1));
bifIndexZlim = [-max(max(abs(bifIndexZlim))) max(max(abs(bifIndexZlim)))];
% plot the topographies
figure
subplot(2,2,1)
cfg = [];
cfg.layout = layoutGTecNew;
cfg.interactive = 'yes';
cfg.zlim = bifIndexZlim*.8;
cfg.xlim = [peakTClust(1,1) peakTClust(2,1)];
cfg.ylim = [peakFreqClust(1,1) peakFreqClust(2,1)];
cfg.marker = 'off';
cfg.comment = 'no';
cfg.interplimits = 'electrodes';
cfg.interpolation = 'v4';
cfg.shading = 'flat';
cfg.gridscale = 250;
ft_topoplotTFR(cfg, dataBifInd)