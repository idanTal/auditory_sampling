%% plot all figures for the auditory sampling paper
clear; close all; clc
%% initialize
subnum = [2 3 4 5 6 7 8 9 10 11 13 14 15 17 18 19 20 21];
subName = {'IT_01', 'ML02', 'ML03', 'MY04', 'ER05', 'CC06', 'LR07', 'IT08', 'SD09', 'AD10', 'SF11', 'NN12', 'CM13', 'EL14','DA15', 'XX16', 'BJ17', 'JS18', 'AS19', 'SY20', 'AG21'};

addpath(genpath('D://Matlab//fieldtrip-20161212'))
addpath(genpath('D://Matlab//CircStat2012a'))
addpath(genpath('D://auditory_sampling//source'))
dataRoot = sprintf('D://auditory_sampling//data');
cd(dataRoot)
% END intialization

%% figure 3 - phase as a function of accuracy


cd(dataRoot);
load phi_clust
% load the indices of the maximal clusters
load coiMax
chns = 1:length(coi);
phiClustAllC = cell(length(chns),length(coi));
phiClustAllE = cell(length(chns),length(coi));
meanPhiClustC = cell(length(coi),1);
meanPhiClustE = cell(length(coi),1);
for cI = 1:length(coi)
    meanPhiClustC{cI} = zeros(length(chns), length(subnum));
    meanPhiClustE{cI} = zeros(length(chns), length(subnum));
    for chnI = 1:length(chns)
        phiClustAllC{chnI,cI} = [];
        phiClustAllE{chnI,cI} = [];
        for subI = 1: length(subnum)
            meanPhiClustC{cI}(chnI,subI) = circ_mean(phiClustC{chnI,subI}(cI,:)');% average phase across trials
            meanPhiClustE{cI}(chnI,subI) = circ_mean(phiClustE{chnI,subI}(cI,:)');% average phase across trials
            tempPhiC = circ_dist(phiClustC{chnI,subI}(cI,:),meanPhiClustC{cI}(chnI,subI));
            tempPhiE = circ_dist(phiClustE{chnI,subI}(cI,:),meanPhiClustC{cI}(chnI,subI));
            phiClustAllC{chnI,cI} = [phiClustAllC{chnI} tempPhiC];
            phiClustAllE{chnI,cI} = [phiClustAllE{chnI} tempPhiE];
        end
    end
end

% 1. bin the trials into 11 bins based on the phase
% 2. identify the trials in each bin and take the average accuracy for
% each subject within a certain bin
% 3. calculate average accuracy across subjects

% start by binning the trials for each subject based on phase (the phase
% information for each subject is in phiClustC (correct) and phiClustE
% (error)

for cI = 1 : length(coi)
    edges = linspace(-pi,pi,12);
    Nratio = zeros(length(subnum), length(chns), length(edges)-1);
    Nc = zeros(length(subnum), length(chns), length(edges)-1);
    Ne = zeros(length(subnum), length(chns), length(edges)-1);
    allPhiC = [];
    allPhiE = [];
    for chnI = 1 : length(chns)
        allPhiC = [];
        allPhiE = [];
        for subI = 1 : length(subnum)
            tempPhiC = circ_dist(phiClustC{chnI,subI}(cI,:),meanPhiClustC{cI}(chnI,subI));
            tempPhiE = circ_dist(phiClustE{chnI,subI}(cI,:),meanPhiClustC{cI}(chnI,subI));
            Nc(subI, chnI, :) = histcounts(tempPhiC, edges);
            Nc(subI, chnI, :) = Nc(subI, chnI, :)./sum(Nc(subI, chnI, :));
            Ne(subI, chnI, :) = histcounts(tempPhiE, edges);
            Ne(subI, chnI, :) = Ne(subI, chnI, :)./sum(Ne(subI, chnI, :));
            Nratio(subI, chnI, :) = (Nc(subI, chnI, :)-Ne(subI, chnI, :))./(Nc(subI, chnI, :)+Ne(subI, chnI, :));
        end
    end
    
    edges = .5*(edges(1:end-1) + edges(2:end));
    meanNratio = squeeze(mean(Nratio,1));
    seNratio = squeeze(std(Nratio,[],1))./sqrt(length(subnum));
    f1 = figure;
    errorbar(edges,nanmean(meanNratio),nanstd(seNratio),'linewidth',1,'color','k')
    set(gca,'xtick',-pi:pi/2:pi)
    set(gca,'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
    set(gca, 'box','off')
    % add a line at 0
    hold on
    line(get(gca,'xlim'),[0 0],'linestyle','--','color',[0.5 0.5 0.5])
    ylabel('(correct-incorrect)/(correct+incorrect)')
    xlabel('Phase bin')

    
    
    %% calculate one way ANOVA with 10 bin levels
    
    % get the average across all channels for the ratio between correct and
    % incorrect responses as a function of the bin
    meanNratioChn = squeeze(nanmean(Nratio,2));
    meanNratioChn(:,ceil(size(meanNratioChn,2)./2)) = [];
    % run ANOVA
    [~,tbl,stats] = anova1(meanNratioChn);
    
    %% check for monotonicity
    h = zeros(1,size(Nratio,3)-1);
    p = zeros(1,size(Nratio,3)-1);
    for binI = 1 : size(Nratio,3)-1
        tempNratio1 = squeeze(nanmean(Nratio(:,:,binI),2));
        tempNratio2 = squeeze(nanmean(Nratio(:,:,binI+1),2));
        diffNratio = tempNratio2-tempNratio1;
        tailSide = 'both';
        [h(binI),p(binI),~,stat] = ttest(tempNratio2,tempNratio1,'Tail', tailSide);
        if stat.tstat < 0 && h(binI) == 1
            h(binI) = -1;
        end
    end
    % correct p-value
    [adjH,~,~,adjP] = fdr_bh(p);
    figure(f1)
    hold on
    tempRatio = nanmean(meanNratio);
    for ii = 1:length(h)
        if adjH(ii) ~= 0
            plot(edges(ii+1)-diff(edges)/2,mean([tempRatio(ii) tempRatio(ii+1)]),'k*')
        end
    end
    title(sprintf('Cluster %d', cI))
end