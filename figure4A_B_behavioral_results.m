%% Analyze behavioral results from auditory sampling experiment
clear; close all; clc
% initialize
dataRoot = 'D://auditory_sampling//data';
subjects = [2,3,4,5,6,7,8,9,10,11,13,14,15,17,18,19,20,21];
% END initialization

respExp = cell(1,length(subjects));
freqExp = cell(1,length(subjects));
accExp = cell(1,length(subjects));
targTimeExp = cell(1,length(subjects));
% load and organize the data without processing
for i = 1 : length(subjects)
    cd(dataRoot)
    % loading behavioral results for the experiment
    load(sprintf('audio_sampling_data%s', num2str(subjects(i)))) 
    
    numTrl = size(data.trlmat_shuffled,2);
    respTemp = data.resp.buttonpress(1,1:numTrl);
    respExp{i} = zeros(1,length(respTemp));
    freqExp{i} = data.trlmat_shuffled(3,1:numTrl);

    % check responses
    for k = 1:length(respTemp)
        if ~iscell(respTemp{k}) && ~isempty(respTemp{k})
            if respTemp{k}(1) == '2' || respTemp{k}(1) == '5'
                respExp{i}(k) = str2double(respTemp{k}(1));
            end
        end
    end
    % check correct trials
    corrTrl = zeros(1,numTrl);
    for k = 1 : numTrl
        if (freqExp{i}(k) == 1 && respExp{i}(k) == 2) || (freqExp{i}(k) == 2 && respExp{i}(k) == 5)
            corrTrl(k) = 1;
        end
    end
    targTimeExp{i} = data.trlmat_shuffled(5,:);
    accExp{i} = corrTrl;
end

%% calculate accuracy for each subject (and SEM)

meanAccExp = zeros(1,length(accExp));
sdAccExp = zeros(1,length(accExp));
for subi = 1:length(accExp)
   meanAccExp(subi) = mean(accExp{subi});
   sdAccExp(subi) = std(accExp{subi})./sqrt(length(accExp{subi}));
end

%% plot error rate as a function of target presentation time
foi = 1:.25:30; % freuqency resolution

% define time bins
tBinSize = 20; % bin size in ms
tBin = 1000:tBinSize:2000; % time of the target stimulus from the onset of the cue

% for each subject, get the ration of correct/incorrect for each bin
targetInd = cell(length(tBin)-1,length(subjects));
accInBin = zeros(length(tBin)-1,length(subjects));
numTrlBin =zeros(length(tBin)-1,length(subjects));
for subI = 1 : length(subjects)
    for tI = 1:length(tBin)-1
        % find indices of targets presented within this bin
        targetInd{tI,subI} = intersect(find(targTimeExp{subI} >= tBin(tI)),find(targTimeExp{subI} < tBin(tI+1)));
        % get the error rate for a single bin
        accInBin(tI,subI) = sum(accExp{subI}(targetInd{tI,subI}) == 0);
        numTrlBin(tI,subI) = length(accExp{subI}(targetInd{tI,subI}));
    end
end
accInBin = accInBin./numTrlBin;

% calculate fft of the binned responses
[x, F] = myPSD(zscore(mean(accInBin,2)'),length(mean(accInBin,2)));

% run permutations
N = 10000;
meanAccBin = mean(accInBin,2);
% correlate a sinosoid with the responses accuracy for each participant
t = 1.02:tBinSize/1000:2;
corrValPerm = zeros(N,round(length(meanAccBin)./2));
for pI = 1 : N
    permI = randperm(length(meanAccBin));
    meanAccBinPerm = meanAccBin(permI);
    [corrValPerm(pI,:), F] = myPSD(zscore(meanAccBinPerm'),length(meanAccBin));
end


% plot the responses and the sin wave
freq2plot = F(maxi(x)); % sin at the frequency of the maximal peak
freq2plotI = find(foi == freq2plot);
figure
% errorbar(t,mean(accInBin,2),std(accInBin,[],2)./sqrt(length(subjects)))
shadedErrorBar(t(1:length(mean(accInBin,2))),smooth(mean(accInBin,2),3),smooth(std(accInBin,[],2),3)./sqrt(length(subjects)),{'color','b', 'linewidth',1})
hold on
tSin = 1.02:.001:2; % time resolution for the sin (for illustration)
sin2plot = ((range(mean(accInBin,2))/3).*sin(2*pi*foi(freq2plotI)*tSin))+mean(mean(accInBin,2));
plot(tSin,sin2plot,'linewidth',1,'color','r')
xlabel('Time from white-noise onset, sec.')
ylabel('error rate')
axis tight
line(get(gca,'xlim'), [mean(mean(accInBin,2)) mean(mean(accInBin,2))])
set(gca,'box','off')
title(sprintf('bin size = %d ms', tBinSize))


% plot the spectrum of the error rate
figure
plot(F,x,'linewidth',1.5)
axis tight
xlabel('Frequency, Hz')
ylabel('Amplitude, a.u.')
xLimit = get(gca,'xlim');
hold on 
prctThr = mean(prctile(abs(corrValPerm),95));
line(xLimit,[prctThr, prctThr] ,'color','k','linewidth',1.5)
set(gca,'box','off')
axis tight
title(sprintf('bin size = %d ms', tBinSize))

