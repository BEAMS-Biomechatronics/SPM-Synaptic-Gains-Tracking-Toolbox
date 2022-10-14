function [RESULT FVAL] = ExhaustedSearch(FeatureVectorsModel,data,windows,config)
% Fit the neural mass model to clinical data by exhausted search. 
%
% [RESULT FVAL] = ExhaustedSearch(FeatureVectorsModel,data,windows,config)
%
% INPUT : 
%       - FeatureVectorsModel: The pre-calculated feature vectors of simulated signals. 
%                              each column corresponds to one triplet indicated  by the first 3 rows
%       - data: EEG data.
%       - windows: output of the function FixWindows, NumOfWin*2, first/second column: start/end of each window in samples.
% 		- config: configuration of algorithm.
%			config.Fs: sample frequency
%           config.NumCandidates: Number of candidates 

% OUTPUT :
%     - RESULT: NumOfWin*1 cells, each cell contains a NumCandidates*3 matrix, which are the first chosen condidates
%     - FVAL: value of error function

RESULT = {};
FVAL = {};
NumWin = size(windows,1);  % number of windows
Params = FeatureVectorsModel(1:3,:);  % all parameters. 
FeatureVectorsModel = FeatureVectorsModel(4:end,:)';   % feature vectors of of all parameter sets, each row is one feature vector
NumCandidates = config.NumCandidates;                  % number of candidates, preset by users. 
try        
    parfor CurWin = 1:NumWin
         datasignal = data(windows(CurWin,1):windows(CurWin,2));          % EEG signal in current window
         N = length(datasignal);
         [FeatureVector]= FeatureCalc(datasignal(:),config);  % feature vector of EEG signal in current window
         %CostAll = pdist2(FeatureVectorsModel,FeatureVector);     % error functions, Euclidean distance between feature vectors
         CostAll = sum((FeatureVectorsModel - repmat(FeatureVector,size(FeatureVectorsModel,1),1)).^2,2);
         [CostAscend, IndexAscend] = sort(CostAll);           % sort the errors
         Candidates = Params(1:3,IndexAscend(1:NumCandidates));   % obtain the candidates
         CostCandidates = CostAll(IndexAscend(1:NumCandidates));  % obtain the corrersponeding error (cost)
         RESULT{CurWin} = Candidates';
         FVAL{CurWin} = CostCandidates';
    end
catch
    for CurWin = 1:NumWin
         datasignal = data(windows(CurWin,1):windows(CurWin,2));          % EEG signal in current window
         N = length(datasignal);
         [FeatureVector]= FeatureCalc(datasignal(:),config);  % feature vector of EEG signal in current window
         CostAll = sum((FeatureVectorsModel - repmat(FeatureVector,size(FeatureVectorsModel,1),1)).^2,2);
         %CostAll = pdist2(FeatureVectorsModel,FeatureVector);     % error functions, Euclidean distance between feature vectors
         [CostAscend, IndexAscend] = sort(CostAll);           % sort the errors
         Candidates = Params(1:3,IndexAscend(1:NumCandidates));   % obtain the candidates
         CostCandidates = CostAll(IndexAscend(1:NumCandidates));  % obtain the corrersponeding error (cost)
         RESULT{CurWin} = Candidates';
         FVAL{CurWin} = CostCandidates';
    end
end
end