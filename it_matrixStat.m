function [pVal, statVal] = it_matrixStat(data, stat)

% This function calculates pixel based statistics on 3D data
% using signrank or ttest against 0

% inputs:
% - data - 3D matrix (e.g. subjects x time x frequency)
% - stat - 'ttest' or 'signrank'

% outputs:
% pVal  - uncorrected p value for each pixel in the matrix
% statVal - summation of the z-value for all pixels in the maximal cluster

%% DEBUGGING
% figure
% % reshape the data for the statistical test
% dataReshape = reshape(data,size(data,1),size(data,2)*size(data,3));
% p = zeros(1,size(dataReshape,2));
% zVal = zeros(1,size(dataReshape,2));
% for ii = 1 : size(dataReshape,2)
%     if ~all(isnan(dataReshape(:,ii)))
%         [p(ii),~,stats] = signrank(dataReshape(:,ii));
%         if ~mod(ii,1000)
%             plot(dataReshape(:,ii))
%             hold on
%         end
%         zVal(ii) = stats.zval;
%     end
% end
% p = reshape(p, size(data,2), size(data,3));
%% END DEBUGGING

% try without reshaping
switch stat
    case 'ttest'
        pVal = nan(size(data,2), size(data,3));
        statVal = nan(size(data,2), size(data,3));
        for ii = 1 : size(data,2)
            for jj = 1 : size(data,3)
                if ~all(isnan(squeeze(data(:,ii,jj))))
                    [~,pVal(ii,jj),~,stats] = ttest(squeeze(data(:,ii,jj)));
                    statVal(ii,jj) = stats.tstat;
                end
            end
        end
    case 'signrank'
        pVal = nan(size(data,2), size(data,3));
        statVal = nan(size(data,2), size(data,3));
        for ii = 1 : size(data,2)
            for jj = 1 : size(data,3)
                if ~all(isnan(squeeze(data(:,ii,jj))))
                    [pVal(ii,jj),~,stats] = signrank(squeeze(data(:,ii,jj)),0);
                    statVal(ii,jj) = stats.zval;
                end
            end
        end
    otherwise
        error('unknown statistical test')
end



