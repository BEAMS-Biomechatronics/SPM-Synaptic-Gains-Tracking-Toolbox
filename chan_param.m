function indAeBG=chan_param(channame,chanlabels)
% look for the channels containing the neuro mass model parameters
[val,tmp,indAeBG]=intersect({[channame '_Ae'],[channame '_B'],[channame '_G']},chanlabels);
if length(indAeBG)~=3
    indAeBG
end