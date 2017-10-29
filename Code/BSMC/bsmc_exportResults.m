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
        load(gtpath); %%load '/home/ethan/Research/Data/ground_truth/cars5_gt.mat';
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
        fprintf(fid,'F-measure based on Avg = %f, AvgMCC = %f\n', fmeasureAvg, mean(MCCV))
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
        load(gtpath); %%load '/home/ethan/Research/Data/ground_truth/cars5_gt.mat';
        for frameNo = start_frame:start_frame+ 9%last_frame
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
        
    case {'tracks'}
        h = vision.TextInserter;
        hm = vision.MarkerInserter;
        colors = [ 0 0 0; 255 0 0; 0 255 0; 0 0 255; 255 255 0; 0 255 255];
        for frameNo = start_frame:last_frame
            % Display motion segmentation output.
            clear mosegState;
            fname = [mosegOutpath '/mosegState-' num2str(frameNo) '.mat'];
            load(fname, 'mosegState');

            lbls = mosegState.lbls;
            outim = imread(sprintf(img_template, frameNo));
            for i=unique(lbls);
                points = round(mosegState.points(1:2, lbls == i));
                hm.Shape = 'Plus';
                hm.BorderColor = 'Custom';
                hm.CustomBorderColor = colors(i+1,:);
                outim = step(hm, outim, int16(points'));
                release(hm);
            end
            
            % Display image number
            h.Text = num2str(frameNo);
            h.FontSize = 24;
            h.Color = [255 255 0]; % Yellow
            h.Location = [50 50];
            outim = step(h, outim);
            
            imwrite(outim,[respath '/' type '-' num2str(frameNo) '.png'],'png');
            release(h);
        end
    case {'seg','bgapp','bgapp_p','fgapp','fgapp_p'}
        for frameNo = start_frame:last_frame
            fname = [outpath '/state-' num2str(frameNo) '.mat'];
            clear state;
            load(fname, 'state');
            if(strcmp(type,'seg'))

                segIm = state.frame;
                segIm = double(segIm); 
                color = [1,0,0; 0,1,0; 0,0,1]; 
                idx = find([state.masks.label] ~= state.bgClust); 
                if(~isempty(idx))
                    cnt = 1; 
                    for fglayer = idx 
                        chnl = find(color(cnt,:)); 
                        segIm(:,:,chnl) = double(segIm(:,:,chnl)) .* ~state.masks(fglayer).mask  + state.masks(fglayer).mask * 255;% *0.4
                        cnt = cnt + 1; 
                    end
            
                end
                segIm = uint8(segIm); 
            else
                switch(type)
                    case 'bgapp'
                        app = bsmc_sampleKDEArray(state.model.bgapp, options);
                    case 'fgapp'
                        app = bsmc_sampleKDEArray(state.model.fgapp, options);
                    case 'bgapp_p'
                        app = bsmc_sampleKDEArray(state.model.bgapp_p, options);
                    case 'fgapp_p'
                        app = bsmc_sampleKDEArray(state.model.fgapp_p, options);
                end
                if strcmp(options.colorspace, 'rgs')
                    segIm = uint8(bsmc_rgs2RGB(app));
                else
                    segIm = uint8(app);
                end
            end
            Text = num2str(frameNo);
            Location = [50 50];        
            outim = insertText(segIm, Location ,Text,'FontSize',24,'BoxColor','white','BoxOpacity',0.4, 'TextColor','yellow');
            imwrite(outim, [respath '/' type '-' num2str(frameNo) '.png'],'png');
        end

   
end
end
