function [angle bbox fidu_XY] = estimatePoseZR(I_Q)
load model3DZhuRamanan Model3D % reference 3D points corresponding to Zhu & Ramanan detections
        angle = 360;
        % detect facial features on query
        addpath(genpath('ZhuRamanan'))
        load('ZhuRamanan/face_p146_small.mat','model');
        model.interval = 5;
        model.thresh = min(-0.65, model.thresh);
            posemap = [];
        if length(model.components)==13 
            posemap = 90:-15:-90;
        elseif length(model.components)==18
            posemap = [90:-15:15 0 0 0 0 0 0 -15:-15:-90];
        else
            error('Can not recognize this model');
            return
        end
        
        if isempty(model)
            return
        end
        I_Q_bs = detect(I_Q, model, model.thresh);
        if isempty(I_Q_bs)
            return
        end
        
        I_Q_bs = clipboxes(I_Q, I_Q_bs);
        I_Q_bs = nms_face(I_Q_bs,0.3);
        bestbs = I_Q_bs(1);
        
        if (isempty(I_Q_bs))
            return;
        end
        x1 = I_Q_bs(1).xy(:,1);
        y1 = I_Q_bs(1).xy(:,2);
        x2 = I_Q_bs(1).xy(:,3);
        y2 = I_Q_bs(1).xy(:,4);
        fidu_XY = [(x1+x2)/2,(y1+y2)/2];
        bbox = [min(fidu_XY(:,1)) min(fidu_XY(:,2)) max(fidu_XY(:,1))-min(fidu_XY(:,1)) max(fidu_XY(:,2))-min(fidu_XY(:,2))];
        angle = posemap(bestbs.c);