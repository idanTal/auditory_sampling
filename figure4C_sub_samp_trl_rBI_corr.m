%% calculate distribution of r values between BI and accuracy when matching the number of trials
clear; close all; clc
%% initialize
subnum = [2 3 4 5 6 7 8 9 10 11 13 14 15 17 18 19 20 21];
subName = {'IT_01', 'ML02', 'ML03', 'MY04', 'ER05', 'CC06', 'LR07', 'IT08', 'SD09', 'AD10', 'SF11', 'NN12', 'CM13', 'EL14','DA15', 'XX16', 'BJ17', 'JS18', 'AS19', 'SY20', 'AG21'};

addpath(genpath('D:\idan\aud_samp_data\source\figureScripts\auditory_sampling\source'))
dataRoot = sprintf('D:\idan\aud_samp_data\source\figureScripts\auditory_sampling\data');

% END intialization

% load the permuted bifurcation index and accuracy
cd(dataRoot)
load BIallSubSampTrials
load meanAccExp
N = size(BIallSubSampTrials,3); % number of subsamplings

% calcualte correlation for each sub-sample of trials
r = zeros(1,N);
p = zeros(1,N);
for permI = 1 : N
    disp(permI)
    meanBItemp = squeeze(mean(BIallSubSampTrials(:,:,permI),2));
    [r(permI), p(permI)] = myCorrPlot(meanAccExp, meanBItemp,'pearson', 0);
end

figure,
histogram(r,50)
axis tight
set(gca,'box','off')
hold on
line([median(r) median(r)], get(gca, 'ylim'), 'linewidth',2, 'color','r')
line([mean(r) mean(r)], get(gca, 'ylim'), 'linewidth',2, 'color','g')
xlabel('Correlation coefficient')
ylabel('Count')


