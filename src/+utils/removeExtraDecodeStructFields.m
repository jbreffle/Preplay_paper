function op = removeExtraDecodeStructFields(resultsStruct_variedDecode, varargin)
% Remove un-necessesary fields from the decoding struct
%
% resultsStruct_variedDecode = utils.removeExtraDecodeStructFields(resultsStruct_variedDecode);


%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'resultsStruct_variedDecode',	@isstruct)
addParameter(inputObj, 'removePmat',  true,  @islogical)
addParameter(inputObj, 'removeRanksVec',  true,  @islogical)
parse(inputObj, resultsStruct_variedDecode, varargin{:});
p = inputObj.Results;


%%
[nParams, nReps] = size(resultsStruct_variedDecode);
ithParam = 1;  ithRep = 1;
nNets = size(resultsStruct_variedDecode(ithParam, ithRep).net, 2);

op = resultsStruct_variedDecode;

for ithParam = 1:nParams
    for ithRep = 1:nReps
        for ithNet = 1:nNets
            try
                % un-needed .replaytrajectory shuffle fields
                op(ithParam, ithRep).net{ithNet}{1}.replaytrajectory = ...
                    rmfield(op(ithParam, ithRep).net{ithNet}{1}.replaytrajectory, ...
                    {'shuffle_rsquare', 'shuffle_maxJump', 'shuffle_slopes', 'shuffle_weightedSlope'});
                if p.removePmat
                    % un-needed .replaytrajectory pMat field
                    op(ithParam, ithRep).net{ithNet}{1}.replaytrajectory = ...
                        rmfield(op(ithParam, ithRep).net{ithNet}{1}.replaytrajectory, ...
                        'pMat');
                end
                if p.removeRanksVec
                    % ranksVec is not needed
                    op(ithParam, ithRep).net{ithNet}{1} = ...
                        rmfield(op(ithParam, ithRep).net{ithNet}{1}, ...
                        'ranksVec');
                end
            end
        end
    end
end

end