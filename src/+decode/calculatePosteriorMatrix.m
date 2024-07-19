function pMat = calculatePosteriorMatrix(trajExpectedSpikes, trajExponExpectedSpikes, spkPerBin)
% pMat = decode.calculatePosteriorMatrix(trajExpectedSpikes, trajExponExpectedSpikes, spkPerBin);
%
% Baye's rules: calculate posterior prob for given traj and event
%


%% Set up

nSpkPerTBin = squeeze(sum(spkPerBin, 3)); % [nTBin x 1] n spikes in tBin


%% Calculate posterior from Bayes

wrking1 = bsxfun(@power, trajExpectedSpikes, spkPerBin); % [nPos x nTbin x nCell]
wrking2 = bsxfun(@rdivide, wrking1, factorial(spkPerBin)); % [nPos x nTbin x nCell]
wrking = bsxfun(@times, wrking2, trajExponExpectedSpikes); % [nPos x nTbin x nCell]
logPost = sum(log(wrking), 3); % log transformed, non-normalized prob [nPos x Tbin]
post = exp(logPost - max(logPost)); % Peak-normalized prob [nPos x Tbin]
post(:, nSpkPerTBin==0) = 0; % so the posterior matrix can be smoothed.
post(isnan(post)) = 0; % remove NaN
pMat = post; % create a posterior matrix for each traj


%% Normalize to columns sum to 1

szPM2 = size(pMat,2);
for i = 1:szPM2 % normalized across positions to 1 for each time bin
    if (sum(pMat(:,i)) > 0)
        pMat(:,i) = pMat(:,i)./sum(pMat(:,i));
    end
end


end