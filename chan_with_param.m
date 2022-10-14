function [chans,IB]=chan_with_param(chanlabels)
tmp={};
for i=1:length(chanlabels)
    tmp{i}=[chanlabels{i} '_Ae'];
end
[val,IA,IB]=intersect(chanlabels,tmp);
chans=chanlabels(IB);