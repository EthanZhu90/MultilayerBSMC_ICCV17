function [ connectTable ] = bsmc_ConnectgtClusters(gtmask, cmask)
    gtlabels = unique(gtmask);
    gtlabels_num = length(gtlabels); 

    clabels = unique(cmask);
    clabels_num = length(clabels); 
    
    costMat = zeros(gtlabels_num, clabels_num); 
    for i = 1:gtlabels_num
        for j = 1:clabels_num
            costMat(i,j) = sum((gtmask(:)==gtlabels(i)) & (cmask(:)== clabels(j))); 
        end
    end
    costMat = -1 * costMat; 
    [assignment, ~] = munkres(costMat);
    connectTable = zeros(gtlabels_num, 2); 
    for i = 1: gtlabels_num
        if(assignment(i) == 0)
            connectTable(i,:) = [gtlabels(i), -1]; %%% too many gt cluster less c cluster.  some gt is not assign a cluster
        else
            connectTable(i,:) = [gtlabels(i), clabels(assignment(i))];
        end
    end
end

