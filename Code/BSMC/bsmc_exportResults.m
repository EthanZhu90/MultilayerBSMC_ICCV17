function [] = bsmc_exportResults(img_template, start_frame, last_frame, mosegOutpath, outpath, respath,gtpath, options, type)

if (~exist(respath,'file'))
    mkdir(respath);
end

switch(type)
        case {'stats2'}
        fid = fopen([respath '/stats2.txt'],'w');
        fprintf(fid, 'Experiment : %s\n', options.description);
        stats = [];
        PrecesionV =[];
        RecallV =[];
        FmeasureV =[];
        MCCV =[];
        Fcount =0;
        load(gtpath); 
        for frameNo = start_frame:last_frame
            
            clear state;
            fname = [outpath '/state-' num2str(frameNo) '.mat'];
            load(fname, 'state');
        
            gim = ground_truth(:,:,frameNo);
            gtmask = (gim > 0);

            fgmask = state.lmask ~= state.bgClust; 
            TP = sum(gtmask(:) & fgmask(:));
            TN = sum(~gtmask(:) & ~fgmask(:));
            FP = sum(~gtmask(:) & fgmask(:));
            FN = sum(gtmask(:) & ~fgmask(:));
            
            precision = TP / (TP + FP);
            recall = TP / (TP + FN);
            if (precision + recall == 0)
                fmeasure = 0;
            else
                fmeasure = 2 * (precision * recall) / (precision + recall);
            end
            
            MCC = (TP*TN - FP*FN) /sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
 
            PrecesionV = [PrecesionV,precision];
            RecallV = [RecallV,recall];
            FmeasureV = [FmeasureV,fmeasure];
            MCCV = [MCCV, MCC];
            Fcount = Fcount+1;
            
            stats = [stats; frameNo TP TN FP FN precision recall fmeasure MCC];
            
            fprintf(fid,'================ Statistics ========================\n');
            fprintf(fid,'Frame %d : TruePositives = %f , TrueNegatives = %f, FalsePositives = %f, FalseNegatives = %f\n', frameNo, TP, TN, FP, FN);
            fprintf(fid,'Frame %d : Precision = %f , Recall = %f, F-measure = %f, MCC = %f\n', frameNo, precision, recall,fmeasure, MCC);
            fprintf(fid,'====================================================\n');
            
        end
        fmeasureAvg = 2 * (mean(PrecesionV) * mean(RecallV)) / (mean(PrecesionV) + mean(RecallV));
        fprintf(fid,'\n\n================ Final Statistics ========================\n');
        fprintf(fid,'Frame Number %d : AvgPrecision = %f , AvgRecall = %f, AvgF-measure = %f\n', Fcount, mean(PrecesionV)...
            , mean(RecallV), mean(FmeasureV));
        fprintf(fid,'F-measure based on Avg = %f, AvgMCC = %f\n', fmeasureAvg, mean(MCCV));
        fprintf(fid,'====================================================\n');
        
        fclose(fid);
        save([respath '/stats2.mat'], 'stats');
        stats = [Fcount, mean(PrecesionV), mean(RecallV), mean(FmeasureV), fmeasureAvg, mean(MCCV)]; 
        save([respath '/stats2_sum.mat'], 'stats');
    
    case {'evalMultiLabel'}
        fid = fopen([respath '/evalMultiLabel.txt'],'w');
        fprintf(fid, 'Experiment : %s\n', options.description);
        stats = [];
        PrecesionV =[];
        RecallV =[];
        FmeasureV =[];
        Fcount =0;
        load(gtpath);
        for frameNo = start_frame:start_frame+ last_frame -1 
            % start_frame:start_frame + 9 if you wanna eval first 10 frames.
            clear state;
            fname = [outpath '/state-' num2str(frameNo) '.mat'];
            load(fname, 'state');
            cmask = state.lmask; 
            gim = ground_truth(:,:,frameNo);
            gtlayer_num = length(unique(gim))-1; 
            connectTable = bsmc_ConnectgtClusters(gim, cmask);
            
            connectTable_tmp = []; 
            for i = 1:size(connectTable, 1)
                if(connectTable(i,1) ~= 0)
                    connectTable_tmp= [connectTable_tmp; connectTable(i,:)];
                end
            end
            connectTable = connectTable_tmp; 
            
            precisionS = zeros(1, gtlayer_num); 
            recallS = zeros(1, gtlayer_num); 
            for i = 1: gtlayer_num
                if(connectTable(i,2) == -1)
                    precisionS(i) = 1;
                    recallS(i) = 0; 
                else
                    precisionS(i) = sum((gim(:) == connectTable(i,1)) & (cmask(:) == connectTable(i,2))) / sum(cmask(:) == connectTable(i,2));
                    recallS(i) = sum((gim(:) == connectTable(i,1)) & (cmask(:) == connectTable(i,2))) / sum(gim(:) == connectTable(i,1)); 
                end
            end
            fmeasureS = 2 * (precisionS.* recallS) / (precisionS + recallS);
            precision = mean(precisionS);
            recall = mean(recallS);
            fmeasure = mean(fmeasureS);
            
            PrecesionV = [PrecesionV, precision];
            RecallV = [RecallV, recall];
            FmeasureV = [FmeasureV, fmeasure];
            Fcount = Fcount+1;
            
            stats = [stats; frameNo precision recall fmeasure];
            
            fprintf(fid,'================ Statistics ========================\n');
            fprintf(fid,'Frame %d : Precision = %f , Recall = %f, F-measure = %f\n', frameNo, precision, recall,fmeasure);
            fprintf(fid,'====================================================\n');
            
        end
        fmeasureAvg = 2 * (mean(PrecesionV) * mean(RecallV)) / (mean(PrecesionV) + mean(RecallV));
        fprintf(fid,'\n\n================ Final Statistics ========================\n');
        fprintf(fid,'Frame Number %d : AvgPrecision = %f , AvgRecall = %f, AvgF-measure = %f\n, F-measure based on Avg = %f\n', Fcount, mean(PrecesionV)...
            , mean(RecallV), mean(FmeasureV), fmeasureAvg);
        fprintf(fid,'====================================================\n');
        
        fclose(fid);
        save([respath '/evalMultiLabel.mat'], 'stats');
        stats = [Fcount, mean(PrecesionV), mean(RecallV), mean(FmeasureV), fmeasureAvg]; 
        save([respath '/evalMultiLabel_sum.mat'], 'stats');
        
   
end
end
