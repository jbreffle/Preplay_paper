function distance = emd(p, q)
% distance = utils.emd(p, q)
%
% Earth mover's distance between probability distributions p and q
% AKA Wasserstein distance
%

%% Checks

% Check if the probability distributions have the same length
if length(p) ~= length(q)
    error('Probability distributions must have the same length.');
end

% Check if the probability distributions sum to 1
%{
if abs(sum(p) - 1) > 1e-6 || abs(sum(q) - 1) > 1e-6
    %error('Probability distributions must sum to 1.');
    distance = nan;
    return
end
%}

%% Calculation

% Calculate the cumulative distribution functions (CDFs)
cdf_p = cumsum(p);
cdf_q = cumsum(q);

% Calculate the EMD
distance = sum(abs(cdf_p - cdf_q));


end
