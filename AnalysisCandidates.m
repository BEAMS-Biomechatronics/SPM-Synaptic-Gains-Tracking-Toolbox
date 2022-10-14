function [RESULT_DBSCAN] = AnalysisCandidates(RESULT,FVAL,windows,config) 
%% Use DBSCAN to obtain one final cluster from the 50 candidate solutions
%Input:
%	- RESULT: 1 * NumWin cells, each is NumOfCandidates * 3 matrix containing the candidates of parameters,
%	          could be the output of ExhaustedSearch.m.
%	- FVAL: 1 * NumWin cells, each is NumOfCandidates * 3 matrix containing the error of candidates,
%	          could be the output of ExhaustedSearch.m.
%	- windows: output of the function FixWindows, NumOfWin*2, first/second column: start/end of each window in samples.
%	- config:configuration of algorithm.
%		config.Fs: sample frequency
%		config.t_step: sengment window step in second
%Output:
%	- RESULT_DBSCAN: struct
%		fields:
%			Param: parameters in the chosen cluster, which has the minimal average error.
%			MeanFval: average error for all window
%			MeanParam: NumWin * 3 matrix, final identified parameters (centroid of the final cluster)
	
% MinPts = 5;
t_step = config.t_step*config.Fs;
N = windows(end,2);
NumWin = length(RESULT);
                for iWin = 1:NumWin
                    Meanfvals = [];
                    fval = FVAL{iWin};
                    CurrentResult = RESULT{iWin};
                    [class] = DBSCAN(RESULT{iWin});
                    NumOfClass = max(class);
                    Result{iWin}.NumOfClass = NumOfClass;
                    for iClass = 1:NumOfClass
                        Index = find(class == iClass);
                        Result{iWin}.class(iClass).NumOfCandidates = length(Index);
                        Result{iWin}.class(iClass).param = CurrentResult(Index,:);
                        Result{iWin}.class(iClass).fval = fval(Index);
                        Result{iWin}.class(iClass).Meanfval = mean(Result{iWin}.class(iClass).fval);
                    end
                    if iClass == 1
                        Result{iWin}.Param = Result{iWin}.class(1).param;
                        Result{iWin}.fval = fval(Index);
                        ClassChosen = 1;
                    else
                        for iClass = 1:NumOfClass
                            Meanfvals(iClass) = Result{iWin}.class(iClass).Meanfval;
                        end
                            ClassChosen = find(Meanfvals==min(Meanfvals));
                            Result{iWin}.Param = Result{iWin}.class(ClassChosen).param;
                            Result{iWin}.fval = Result{iWin}.class(ClassChosen).fval;
                    end
                    %% average Ae, B and G, instead of average the ratios
                    RESULT_DBSCAN.Param{iWin} = Result{iWin}.Param;
                    RESULT_DBSCAN.MeanFval(iWin) = mean(Result{iWin}.fval);
                    Ae_DBSCAN{iWin} = RESULT_DBSCAN.Param{iWin}(:,1);
                    B_DBSCAN{iWin} = RESULT_DBSCAN.Param{iWin}(:,2);
                    G_DBSCAN{iWin} = RESULT_DBSCAN.Param{iWin}(:,3);
                %% indicators
                    MeanAe_DBSCAN(iWin) = mean(Ae_DBSCAN{iWin});
                    MeanB_DBSCAN(iWin) = mean(B_DBSCAN{iWin});
                    MeanG_DBSCAN(iWin) = mean(G_DBSCAN{iWin});
                   
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                RESULT_DBSCAN.MeanParam = [MeanAe_DBSCAN;MeanB_DBSCAN; MeanG_DBSCAN]'; %%final mean of Ae B G;

end


