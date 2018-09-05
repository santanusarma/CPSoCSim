function Labels=RemoveUnderscoreLabels(trace)
%Reformat the underscores in the labes for legend in plot

for i=1:length(trace.Labels)
    Labels(i)=strrep(trace.Labels(i), '_', '\_')
end