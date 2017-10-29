function [] = sumResult(  )
segname_vector= {'cars1','cars2','cars3','cars4','cars5','cars6','cars7','cars8','people1','people2'}; 
stats_store = []; 
for i = 1: length(segname_vector)
    outpath = sprintf('../Results/Results/%s/', segname_vector{i});
    statsfile  = [outpath, 'stats2_sum.mat'];  % 'evalMultiLabel_sum.mat']; 
    load(statsfile); 
    stats_store = [stats_store;  stats]; 
end
stats_store
stats_Avg = mean(stats_store, 1); 

disp 'The average of stats:'
stats_Avg

frame_Num = sum(stats_store(:,1)); 

stats_tmp = zeros(1,4); 
for i = 1: length(segname_vector)
   stats_tmp = stats_tmp+ stats_store(i,2:5).* stats_store(i,1); 
end
stats_tmp = stats_tmp./frame_Num; 
stats_Avg(2:5) =  stats_tmp; 

disp 'The average of stats with frame weighted:'
stats_Avg

end

