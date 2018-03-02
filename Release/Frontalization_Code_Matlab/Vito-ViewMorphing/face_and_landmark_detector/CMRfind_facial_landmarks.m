function landmarks = CMRfind_facial_landmarks(X,bbox, Ldetector,scfak)

% Init

[a,b,c] = size(X);
if c~=1
   X=uint8(mean(X,3)); 
%    X=uint8(X(:,:,1)); 
end
nfaces = size(bbox,1);
bbox1 = bbox;
landmarks = zeros(nfaces,136);

for k1=1:nfaces;
    bbox = bbox1(k1,:);
    
    % find facial region and rescale image
    new_scale_fak = Ldetector.bbox_size/bbox(3);
    R=imresize(X,[round(a*new_scale_fak), round(b*new_scale_fak)],'bilinear');
    bbox = round(bbox*new_scale_fak);

    % find initial shape position
    ref_position.y1 = round(bbox(2)+Ldetector.init_position(1)*bbox(4));
    ref_position.y2 = round(bbox(2)+scfak*bbox(4));
    ref_position.x1 = round(bbox(1)+Ldetector.init_position(2)*bbox(3)); 
    ref_position.x2 = round(bbox(1)+Ldetector.init_position(3)*bbox(3)); 
    x_diff = abs(ref_position.x1-ref_position.x2);
    y_diff = abs(ref_position.y1-ref_position.y2);

    %get means coordinates
    live_position.y1 = round((Ldetector.init_model(38,2)+Ldetector.init_model(45,2))/2); %points 38 and 45 represent the eyes
    live_position.y2 = round(Ldetector.init_model(9,2)); %point 9 represents the chin
    live_position.x1 = round(Ldetector.init_model(38,1)); 
    live_position.x2 = round(Ldetector.init_model(45,1));
    live_scale_fak_x = abs(live_position.x1-live_position.x2)/x_diff;
    live_scale_fak_y = abs(live_position.y1-live_position.y2)/y_diff;

    ref_point1 = [live_position.x1/live_scale_fak_x, live_position.y1/live_scale_fak_y];
    displacement = [ref_position.x1,ref_position.y1] - ref_point1;
    avg_adjusted = zeros(68,2);
    for j=1:68
        avg_adjusted(j,:) = [Ldetector.init_model(j,1)/live_scale_fak_x+displacement(1,1),Ldetector.init_model(j,2)/live_scale_fak_y+displacement(1,2)];
    end
    init_position_vect = avg_adjusted(:);


    % find landmarks 
    interest_points = reshape(round(init_position_vect),[68,2]);
    for j=1:Ldetector.niter 
        if scfak==1
            [features1, noValidPoints] = MyextractHOGFeatures1(R, round(interest_points), Ldetector.csiz(j), Ldetector.bsiz(j), 0, 9, 0);
        else
            [features1, noValidPoints] = MyextractHOGFeatures1(R, round(interest_points), Ldetector.csiz(j), Ldetector.bsiz(j), 0, 9, 1);
        end
       if noValidPoints~=68
           R_new = padarray(R,[double(round(2.5*32)) double(round(2.5*32))],'replicate');
           interest_points1 = reshape((interest_points(:))+(2.5*32),[136/2,2]);
           if scfak==1
               [features1, noValidPoints] = MyextractHOGFeatures1(R_new, round(interest_points1), Ldetector.csiz(j), Ldetector.bsiz(j), 0, 9, 0);   
           else
               [features1, noValidPoints] = MyextractHOGFeatures1(R_new, round(interest_points1), Ldetector.csiz(j), Ldetector.bsiz(j), 0, 9, 0);  
           end
               
       end
        
        if Ldetector.bsiz(j)==4
            features_t = double(Ldetector.model4.W)'*(features1(:)-double(Ldetector.model4.P));
        elseif Ldetector.bsiz(j)==2
            features_t = double(Ldetector.model2.W)'*(features1(:)-double(Ldetector.model2.P));
        end
       sims = evaluate_gmm1(double(features_t),double(Ldetector.gmm(j).Mean),double(Ldetector.gmm(j).Cov));
       Xf=[features_t;1];
       tmp=(sims(1))*double(Ldetector.gmm(j).transform(1).W');
       for k=2:size(Ldetector.gmm(j).transform,2)      
           tmp=tmp+(sims(k))*double(Ldetector.gmm(j).transform(k).W');
       end
       tmp = tmp*Xf;
       
       new_positions = (interest_points(:) + tmp);
       interest_points = reshape(new_positions,[136/2,2]);
    end 
    tmp = round(interest_points/new_scale_fak);
    landmarks(k1,:) = tmp(:);
end

end