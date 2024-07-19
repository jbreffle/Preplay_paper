function gridDispProgress(~, totalCount)
% gridDispProgress.m
%
% Prints the progress and estimated remaining time, even in a parfor loop
%
% Example Usage:
%
% D = parallel.pool.DataQueue;
% num_files = 10;
% gridDispProgress(1, num_files);
% afterEach(D, @gridDispProgress);
% tic
% parfor ithParamSet = 1:num_files
%     send(D, 1);
% end
%

persistent TOTAL COUNT
if nargin == 2
    % initialisation mode
    TOTAL = totalCount;
    COUNT = 0;
    disp(['Starting parameter grid with ', num2str(TOTAL), ' parameter points'])
else
    % afterEach call, increment COUNT
    COUNT = 1 + COUNT;
    fracCompleted = COUNT / TOTAL;
    disp(['Finished parameter set ', num2str(COUNT), ' of ', num2str(TOTAL), ...
        newline, '    Total runtime: ', num2str(toc/60/60), ' hours', ...
        newline, '    Est. remaining time: ', num2str(( (toc/fracCompleted)-toc )/60/60), ' hours'
        ])
end
end