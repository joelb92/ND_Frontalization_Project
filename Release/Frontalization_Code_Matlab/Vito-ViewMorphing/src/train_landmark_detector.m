% SDM with pre-initialization

clc
clear
%close all

%addpath('I:\detection\300w\eng\sup')
addpath('C:\opencv\mexopencv')
extractor = cv.DescriptorExtractor('SURF'); 

% Parameters
permuts = 5; %location permutations of the initial points
niter = 6; %number of iterations
kscale = linspace(1/10,1/20,niter);%keypoint scale (HOG: 1/18, SURF: 1/6)
blockSize = 4;
lambda = linspace(10,1,niter); 
verb = 1; 

%% Parse training data

data_dir = {'I:\detection\300w\annotations\lfpw\trainset\', 'I:\detection\300w\annotations\ibug\', 'I:\detection\300w\annotations\afw\', 'I:\detection\300w\annotations\helen\trainset\'};
points_list = {cellstr(ls([data_dir{1},'*.pts'])), cellstr(ls([data_dir{2},'*.pts'])), cellstr(ls([data_dir{3},'*.pts'])), cellstr(ls([data_dir{4},'*.pts']))};
image_list =  {cellstr(ls([data_dir{1},'*.png'])), cellstr(ls([data_dir{2},'*.jpg'])), cellstr(ls([data_dir{3},'*.jpg'])), cellstr(ls([data_dir{4},'*.jpg']))};

t = 1;
for d = 1 : length(data_dir)
    for c = 1 : length(points_list{d})
        train(t).path = data_dir{d};
        train(t).name = image_list{d}{c};
        train(t).fp = ptsread([data_dir{d},points_list{d}{c}])';
        t = t + 1;
    end
end

%% Points of an average face

data_avg = mean(cat(3, train.fp),3);
data_avg = data_avg - repmat(mean(data_avg,2),1,size(data_avg,2));

% Get point coordinates from the average face
live_position.y1 = round(mean(data_avg(2,37:48))); %points 37 to 48 represent the eyes
live_position.y2 = round(data_avg(2,9)); %point 9 represents the chin
live_position.x1 = round(mean(data_avg(1,37:42))); %right eye
live_position.x2 = round(mean(data_avg(1,43:48))); %left eye

if verb
    fh = figure;
    lbls = cellstr(num2str([1:size(data_avg,2)]'));
    plot(data_avg(1,:),data_avg(2,:),'.')
    hold on
    set(gca,'YDir','reverse')
    text(data_avg(1,:), data_avg(2,:), lbls, 'VerticalAlignment','bottom', ...
                             'HorizontalAlignment','right')
    plot(live_position.x1,live_position.y1,'ro','MarkerFaceColor','r','MarkerSize',8)
    plot(live_position.x2,live_position.y1,'go','MarkerFaceColor','g','MarkerSize',8)
    plot(live_position.x2,live_position.y2,'bo','MarkerFaceColor','b','MarkerSize',8)
    set(fh,'Name','Fiducial points of an average face')
end

%% Face detection (Viola-Jones detector) and setting of the initial locations on training images

y_poss = 0.39;  %Delez ipsilona v odstotkih glede na velikost kvadrata obraza za oci
x_pos1 = 0.34;  %Delez iksa v odstotkih glede na velikost kvadrata obraza za levo oko
x_pos2 = 0.66;  %Delez iksa v odstotkih glede na velikost kvadrata obraza za desno oko

mdu = []; %images where no faces were detected

if verb
    fh = figure;
    set(fh, 'Name', 'Face detection on training images')
end

for c = 1 : length(train)
    
    img = imread([train(c).path,train(c).name]);
    if ~isa(img,'uint8')
        img = im2uint8(img);
    end
    rect = face_detector(img);
    
    if ~isempty(rect)
        % Compute reference coordinates - aproximate eye and chin locations
        % based on the face rectangle location
        ref_position.x1 = round(rect(1)+x_pos1*rect(3)); 
        ref_position.y1 = round(rect(2)+y_poss*rect(4)); 
        ref_position.x2 = round(rect(1)+x_pos2*rect(3)); 
        ref_position.y2 = round(rect(2)+rect(4));
        
        x_diff = abs(ref_position.x1-ref_position.x2);
        y_diff = abs(ref_position.y1-ref_position.y2);
        live_scale_fak_x = abs(live_position.x1-live_position.x2)/x_diff;
        live_scale_fak_y = abs(live_position.y1-live_position.y2)/y_diff;
        live_scale_fak = mean([live_scale_fak_x,live_scale_fak_y]);
        ref_point1 = [live_position.x1/live_scale_fak, live_position.y1/live_scale_fak];
        
        % Compute displacement and scaling
        displacement = [ref_position.x1,ref_position.y1] - ref_point1;
        scale_var = 0.05;
        disp_var = 0.05*rect(3);
        samples = normrnd(repmat([live_scale_fak displacement(1) displacement(2)],permuts,1),...
                repmat([scale_var*live_scale_fak disp_var disp_var],permuts,1),[permuts,3]);

        % Recompute points
        for k = 1 : permuts
            avg_adjusted = zeros(size(data_avg));
            for j = 1 : size(avg_adjusted,2)
                avg_adjusted(:,j) = [data_avg(1,j)/samples(k,1)+samples(k,2),...
                                     data_avg(2,j)/samples(k,1)+samples(k,3)];
            end
            train(c).fp0{k} = avg_adjusted;
        end

        % Remeber results
        train(c).rect = [rect(1),rect(2),rect(3),rect(4)]; % koordinate kvadrata obraza [x,y,w,h]
        
        fprintf('Done with face detection in a training image "%s" (#%d/%d).\n',train(c).name,c,length(train));

        if verb
            hold off
            imshow(img,[]); hold on
            rectangle('Position',[rect(1),rect(2),rect(3),rect(4)],'EdgeColor','red')
            plot(ref_position.x1,ref_position.y1,'ro','MarkerFaceColor','r','MarkerSize',8)
            plot(ref_position.x2,ref_position.y1,'go','MarkerFaceColor','g','MarkerSize',8)
            plot(ref_position.x2,ref_position.y2,'bo','MarkerFaceColor','b','MarkerSize',8)
            for k = 1 : permuts
                plot(train(c).fp0{k}(1,:),train(c).fp0{k}(2,:),'.');
            end
            plot(train(c).fp(1,:),train(c).fp(2,:),'y+')
            drawnow; %pause
        end
    else
        fprintf('Unable to detect any face on a training image "%s" (#%d/%d).\n',train(c).name,c,length(train))
        mdu = [mdu;c];
    end
end
train(mdu)=[];

%% Find misdetected training images 
% Criterion: more than half of hand labeled points are in the rectangle area

fp_err = [];
for c = 1 : length(train)
    if sum(train(c).rect(1) < train(c).fp(1,:) &...
        train(c).rect(2) < train(c).fp(2,:)&...
       (train(c).rect(1) + train(c).rect(3)) > train(c).fp(1,:)&...
       (train(c).rect(2) + train(c).rect(4)) > train(c).fp(2,:)) < length(train(c).fp)*.5 || ...
        prod(max(train(c).fp,[],2)-min(train(c).fp,[],2))/prod(max(train(c).fp0{1},[],2)-min(train(c).fp0{1},[],2)) > 3 || ...
        prod(max(train(c).fp,[],2)-min(train(c).fp,[],2))/prod(max(train(c).fp0{1},[],2)-min(train(c).fp0{1},[],2)) < 1/3
            fp_err = [fp_err c];
    end
end

if verb
    fh = figure;
    set(fh,'Name','Misdetected training faces.')
    for c = fp_err
        img = imread([train(c).path,train(c).name]);
        imshow(img,[])
        hold on
        rectangle('Position',[train(c).rect(1),train(c).rect(2),train(c).rect(3),train(c).rect(4)],'EdgeColor','red')
        plot(train(c).fp(1,:),train(c).fp(2,:),'y+')
        hold off
        drawnow
        %pause
    end
end
train(fp_err) = [];
fprintf('\n Removed %d misdetected training images.\n',length(fp_err))

% =========================================================================
%% Train descent directions

[train.fpt] = train.fp0;

if verb
    hf = figure('units','normalized','outerposition',[0 0 1 1]);
    set(hf, 'Name', 'Point convergence on a randomly selected training image')
    n = randi(length(train));
end

y_ref = zeros(permuts*length(train),size(train(1).fp,2)*2);
for c = 1 : length(train)
    for p = 1 : permuts
        y_ref((c-1)*permuts+p,:) = train(c).fp(:);
    end
end

for k = 1 : niter

    y_pos = zeros(permuts*length(train),size(train(1).fp,2)*2);
    for c = 1 : length(train)
        for p = 1 : permuts
            y_pos((c-1)*permuts+p,:) = reshape(train(c).fpt{p},[],1);
        end
    end

    fkt = []; fks = [];

    % Descriptors
    for c = 1:length(train)

        img = imread([train(c).path,train(c).name]);
        if ~isa(img,'uint8')
            img = im2uint8(img);
        end
        img = uint8(normalize8(mean(img,3)));
        
        for p = 1 : permuts

            % HOG (Matlab)
            kpt_scale = round(kscale(k)*mean([train(c).rect(3) train(c).rect(4)]));
            %[fks{(c-1)*permuts+p},~] = extractHOGFeatures(img,train(c).fpt{p}','CellSize',[kpt_scale kpt_scale],'BlockSize',[blockSize blockSize]);
            cvkpt = cell2struct([num2cell(train(c).fpt{p}',2), repmat(num2cell(kpt_scale),size(train(c).fpt{p}',1),1)], [{'pt'};{'size'}], 2); 
            fks{(c-1)*permuts+p} = extractor.compute(img, cvkpt);
        end

    end

    rema = cellfun('size',fks,1)==length(train(1).fp);
    fks = cellfun(@(x)x(:),fks(rema),'UniformOutput',false);
    fks = cell2mat(fks); fks = double(fks);
    %[coeff,score,latent,tsquared,explained,mu] = pca(fks');
    %pcn = cumsum(explained) < 98;
    %coeff = coeff(:,pcn);
    fk_mean_s{k} = mean(fks,2);
    fks = fks - repmat(fk_mean_s{k},1,size(fks,2));
    %fks = coeff*(fks'*coeff)';
    fks = fks';
    xs = y_ref(rema,:);
    xk = y_pos(rema,:);
    dxk = (xs - xk);
    dxk_mean{k} = mean(dxk,1);
    dxk = dxk - repmat(dxk_mean{k},size(dxk,1),1);

    if 0 %ridge regression
        Rb{k} = (fks'*fks + lambda(k)*eye(size(fks,2)))\fks' * dxk;
        xk = xk + fks*Rb{k} + repmat(dxk_mean{k},size(dxk,1),1);
    else %kernel ridge regression
        K{k} = kernelmatrix('poly',fks',[],0,3);
        X{k} = fks; Y{k} = dxk;
        xk = xk + [Y{k}'*((K{k}+lambda(k)*eye(size(Y{k},1)))\K{k})]' + repmat(dxk_mean{k},size(Y{k},1),1);
    end


    % Plot
    if verb
        fh1 = subplot(2,ceil(niter/2),k);
        img = imread([train(n).path,train(n).name]);
        fh2 = imshow(img); hold on
        fptp = train(n).fpt{1};
        axis equal
        plot(train(n).fpt{1}(1,:),train(n).fpt{1}(2,:),'r.')
        plot(train(n).fp(1,:),train(n).fp(2,:),'g.');
        title(fh1,sprintf('%d. iteration',k))
        drawnow
    end
    fprintf('\tIteration #%d.\n',k)

    t1 = 1; t2 = 1;
    for c = 1 : length(rema)
        if rema(c)
            train(ceil(c/permuts)).fpt{t2} = reshape(xk(t1,:),2,size(train(1).fp,2));
            t1 = t1 + 1;
        end
        t2 = t2 + 1;
        if t2 > permuts, t2 = 1; end
    end

end

%==========================================================================
%% Train pre-initializer
% 
% % Use AFW (data set)
% pini_dir = 'I:\detection\300w\annotations\ibug\';
% pts_list = cellstr(ls([pini_dir,'*.pts']));
% figs_list = cellstr(ls([pini_dir,'*.jpg']));
% 
% 
% % Parse points
% pini = [];
% for c = 1 : length(pts_list)
%     pini(c).fp = ptsread([pini_dir,pts_list{c}])';
%     pini(c).name = figs_list{c};
% end
% 
% %% Face detection (Viola-Jones detector) and setting of the initial locations on preinit images
% 
% mdu = []; %images where no faces were detected
% 
% if verb
%     fh = figure;
%     set(fh, 'Name', 'Face detection on pre-init images')
% end
% 
% %pini = rmfield(pini,'fp0');
% 
% for c = 1 : length(pini)
%     
%     img = imread([pini_dir,pini(c).name]);
%     if ~isa(img,'uint8')
%         img = im2uint8(img);
%     end
%     rect = face_detector(img);
%     
%     if ~isempty(rect)
%         % Compute reference coordinates - aproximate eye and chin locations
%         % based on the face rectangle location
%         ref_position.x1 = round(rect(1)+x_pos1*rect(3)); 
%         ref_position.y1 = round(rect(2)+y_poss*rect(4)); 
%         ref_position.x2 = round(rect(1)+x_pos2*rect(3)); 
%         ref_position.y2 = round(rect(2)+rect(4));
%         
%         x_diff = abs(ref_position.x1-ref_position.x2);
%         y_diff = abs(ref_position.y1-ref_position.y2);
%         live_scale_fak_x = abs(live_position.x1-live_position.x2)/x_diff;
%         live_scale_fak_y = abs(live_position.y1-live_position.y2)/y_diff;
%         live_scale_fak = mean([live_scale_fak_x,live_scale_fak_y]);
%         ref_point1 = [live_position.x1/live_scale_fak, live_position.y1/live_scale_fak];
%         displacement = [ref_position.x1,ref_position.y1] - ref_point1;
%         
%         % Recompute points
%         scl = normrnd(1,.2,1,15);
%         dsplc = normrnd(0,.1,2,15);
%         for t = 1 : 15
%             pini(c).fp0{t} = data_avg./live_scale_fak.*scl(t) +...
%                           repmat(dsplc(:,t).*rect(3:4)',1,size(data_avg,2)) +...
%                           repmat(displacement',1,size(data_avg,2));
%             pini(c).sd(:,t) = [scl(t); dsplc(:,t)];
%             t = t + 1;
%         end
%         
%         % Remeber results
%         pini(c).rect = [rect(1),rect(2),rect(3),rect(4)];
%         
%         fprintf('Done with face detection in a pre-init image "%s" (#%d/%d).\n',pini(c).name,c,length(pini));
% 
%         if verb
%             hold off
%             imshow(img,[]); hold on
%             rectangle('Position',[rect(1),rect(2),rect(3),rect(4)],'EdgeColor','red')
%             plot(ref_position.x1,ref_position.y1,'ro','MarkerFaceColor','r','MarkerSize',8)
%             plot(ref_position.x2,ref_position.y1,'go','MarkerFaceColor','g','MarkerSize',8)
%             plot(ref_position.x2,ref_position.y2,'bo','MarkerFaceColor','b','MarkerSize',8)
%             for t = 1 : length(pini(c).fp0)
%                 plot(pini(c).fp0{t}(1,:),pini(c).fp0{t}(2,:),'.');
%             end
%             plot(pini(c).fp(1,:),pini(c).fp(2,:),'y+')
%             drawnow; %pause
%         end
%     else
%         fprintf('Unable to detect any face on a pre-init image "%s" (#%d/%d).\n',pini(c).name,c,length(pini))
%         mdu = [mdu;c];
%     end
% end
% pini(mdu)=[];
% 
% %% Find misdetected (pre-initialization) images 
% % Criterion: more than a half of the hand labeled points are in the rectangle area
% 
% fp_err = [];
% for c = 1 : length(pini)
%     if sum(pini(c).rect(1) < pini(c).fp(1,:) &...
%         pini(c).rect(2) < pini(c).fp(2,:)&...
%        (pini(c).rect(1) + pini(c).rect(3)) > pini(c).fp(1,:)&...
%        (pini(c).rect(2) + pini(c).rect(4)) > pini(c).fp(2,:)) < length(pini(c).fp)*.5 || ...
%         prod(max(pini(c).fp,[],2)-min(pini(c).fp,[],2))/prod(max(pini(c).fp0{1},[],2)-min(pini(c).fp0{1},[],2)) > 3 || ...
%         prod(max(pini(c).fp,[],2)-min(pini(c).fp,[],2))/prod(max(pini(c).fp0{1},[],2)-min(pini(c).fp0{1},[],2)) < 1/3;
%             fp_err = [fp_err c];
%     end
% end
% 
% if verb && ~isempty(fp_err)
%     fh = figure;
%     set(fh,'Name','Misdetected pre-init faces.')
%     for c = fp_err
%         img = imread([pini_dir,pini(c).name]);
%         imshow(img,[])
%         hold on
%         rectangle('Position',[pini(c).rect(1),pini(c).rect(2),pini(c).rect(3),pini(c).rect(4)],'EdgeColor','red')
%         plot(pini(c).fp(1,:),pini(c).fp(2,:),'y+')
%         hold off
%         drawnow
%         %pause
%     end
% end
% pini(fp_err) = [];
% fprintf('\n%d misdetected pre-init images.\n',length(fp_err))
% 
% %% SDM on pre-init images
% 
% [pini.fpt] = pini.fp0;
% 
% t = 1;
% clear tsr M fks_pi param_pi fks fks_pi fks_pi2
% 
% for c = 1 : length(pini)
% 
%     img = imread([pini_dir,pini(c).name]);
%     if ~isa(img,'uint8')
%         img = im2uint8(img);
%     end
%     img = uint8(normalize8(mean(img,3)));
%     
%     [xi,yi] = meshgrid(linspace(pini(c).rect(1),pini(c).rect(1)+pini(c).rect(3),10),...
%                        linspace(pini(c).rect(2),pini(c).rect(2)+pini(c).rect(4),10));
%     pts = [xi(:),yi(:)];
%     
%     kpt_scale_pi =  round(kscale(1)*mean([pini(c).rect(3) pini(c).rect(4)]));
% %     [fks_pi,~] = extractHOGFeatures(img, pts,'CellSize',[kpt_scale kpt_scale],'BlockSize',[blockSize blockSize]);
%     cvkpt = cell2struct([num2cell(pts,2), repmat(num2cell(kpt_scale_pi),size(pts,1),1)], [{'pt'};{'size'}], 2); 
%     fks_pi = extractor.compute(img, cvkpt);
%     
%     if size(fks_pi,1) < length(pts)
%         border_size = 100;
%         img_tmp = [zeros(size(img,1),border_size),img,zeros(size(img,1),border_size)];
%         img_tmp = [zeros(border_size,size(img_tmp,2));img_tmp;zeros(border_size,size(img_tmp,2))];
%         %[fks_pi,~] = extractHOGFeatures(img_tmp, pts+border_size,'CellSize',[kpt_scale kpt_scale],'BlockSize',[blockSize blockSize]);
%         cvkpt = cell2struct([num2cell(pts+border_size,2), repmat(num2cell(kpt_scale_pi),size(pts,1),1)], [{'pt'};{'size'}], 2); 
%         fks_pi = extractor.compute(img_tmp, cvkpt);
%     end
%     
%     pts_err = [];
%     for p = 1 : length(pini(c).fp0)
%         
%         fpt = pini(c).fp0{p};
%         
%         for k = 1 : niter
%             
%             kpt_scale =  round(kscale(k)*mean([pini(c).rect(3) pini(c).rect(4)]));
% %             [fks,~] = extractHOGFeatures(img, fpt','CellSize',[kpt_scale kpt_scale],'BlockSize',[blockSize blockSize]);
%             cvkpt = cell2struct([num2cell(fpt',2), repmat(num2cell(kpt_scale),size(fpt',1),1)], [{'pt'};{'size'}], 2); 
%             fks = extractor.compute(img, cvkpt);
%             if size(fks,1) < size(data_avg,2)
%                 border_size = 100;
%                 img_tmp = [zeros(size(img,1),border_size),img,zeros(size(img,1),border_size)];
%                 img_tmp = [zeros(border_size,size(img_tmp,2));img_tmp;zeros(border_size,size(img_tmp,2))];
%                 %[fks,~] = extractHOGFeatures(img_tmp, fpt'+border_size,'CellSize',[kpt_scale kpt_scale],'BlockSize',[blockSize blockSize]);
%                 cvkpt = cell2struct([num2cell(fpt'+border_size,2), repmat(num2cell(kpt_scale),size(fpt',1),1)], [{'pt'};{'size'}], 2); 
%                 fks = extractor.compute(img_tmp, cvkpt);
%             end
%             fks = (fks(:) - fk_mean_s{k})';
% 
%             if 1
%                 tmp = fpt(:)' + fks*Rb{k} + dxk_mean{k};
%             else
%                 tmp = fpt(:)' + (Y{k}'*((K{k}+lambda(k)*eye(size(Y{k},1)))\kernelmatrix('poly',X{k}',fks',0,3)))' + dxk_mean{k};
%             end
% 
%             fpt = reshape(tmp,2,[]);
%         end
%         pts_err(p) = sqrt(sum((pini(c).fp(:) - fpt(:)).^2))/length(pini(c).fp);
%     end
%     
%     idx = find(pts_err == min(pts_err));
%     fks_pi2(c,:) = fks_pi(:); 
%     param_pi(c,:) = pini(c).sd(:,idx);
%     disp(c)
% end
% fks_pi2_mean = mean(fks_pi2);
% fks_pi3 = fks_pi2 - repmat(fks_pi2_mean,size(fks_pi2,1),1);
% param_pi_mean = mean(param_pi);
% param_pi1 = param_pi - repmat(param_pi_mean,size(param_pi,1),1);
% R_pi = (fks_pi3'*fks_pi3 + 5*eye(size(fks_pi3,2)))\fks_pi3' * param_pi1;

%==========================================================================
%% Parse test data

data_dir = 'I:\detection\300w\annotations\afw\';
points_list = cellstr(ls([data_dir,'*.pts']));
image_list = cellstr(ls([data_dir,'*.jpg']));

clear test
for c = 1 : length(points_list)
    test(c).path = data_dir;
    test(c).name = image_list{c};
    test(c).fp = ptsread([data_dir,points_list{c}])';
end

%% Face detection (Viola-Jones detector) and setting of the initial locations on test images

mdu = []; %images where no faces were detected

if verb
    fh = figure;
    set(fh, 'Name', 'Face detection on test images')
end

for c = 1 : length(test)
    
    img = imread([test(c).path,test(c).name]);
    if ~isa(img,'uint8')
        img = im2uint8(img);
    end
    rect = face_detector(img);
    
    if ~isempty(rect)
        % Compute reference coordinates - aproximate eye and chin locations
        % based on the face rectangle location
        ref_position.x1 = round(rect(1)+x_pos1*rect(3)); 
        ref_position.y1 = round(rect(2)+y_poss*rect(4)); 
        ref_position.x2 = round(rect(1)+x_pos2*rect(3)); 
        ref_position.y2 = round(rect(2)+rect(4));
        
        x_diff = abs(ref_position.x1-ref_position.x2);
        y_diff = abs(ref_position.y1-ref_position.y2);
        live_scale_fak_x = abs(live_position.x1-live_position.x2)/x_diff;
        live_scale_fak_y = abs(live_position.y1-live_position.y2)/y_diff;
        live_scale_fak = mean([live_scale_fak_x,live_scale_fak_y]);
        ref_point1 = [live_position.x1/live_scale_fak, live_position.y1/live_scale_fak];
        displacement = [ref_position.x1,ref_position.y1] - ref_point1;

        % Recompute points
        for j = 1 : size(data_avg,2)
            test(c).fp0(:,j) = [data_avg(1,j)/live_scale_fak + displacement(1),...
                                data_avg(2,j)/live_scale_fak + displacement(2)];
        end


        % Remeber results
        test(c).rect = [rect(1),rect(2),rect(3),rect(4)];
        
        fprintf('Done with face detection in a test image "%s" (#%d/%d).\n',test(c).name,c,length(test));

        if verb
            hold off
            imshow(img,[]); hold on
            rectangle('Position',[rect(1),rect(2),rect(3),rect(4)],'EdgeColor','red')
            plot(ref_position.x1,ref_position.y1,'ro','MarkerFaceColor','r','MarkerSize',8)
            plot(ref_position.x2,ref_position.y1,'go','MarkerFaceColor','g','MarkerSize',8)
            plot(ref_position.x2,ref_position.y2,'bo','MarkerFaceColor','b','MarkerSize',8)
            plot(test(c).fp0(1,:),test(c).fp0(2,:),'r.');
            plot(test(c).fp(1,:),test(c).fp(2,:),'g.')
            drawnow; %pause
        end
    else
        fprintf('Unable to detect any face on a test image "%s" (#%d/%d).\n',test(c).name,c,length(test))
        mdu = [mdu;c];
    end
end
test(mdu)=[];

%% Find misdetected test images 
% Criterion: more than 50% of hand labeled points are in the rectangle area

fp_err = [];
for c = 1 : length(test)
    if sum(test(c).rect(1) < test(c).fp(1,:) &...
        test(c).rect(2) < test(c).fp(2,:)&...
       (test(c).rect(1) + test(c).rect(3)) > test(c).fp(1,:)&...
       (test(c).rect(2) + test(c).rect(4)) > test(c).fp(2,:)) < length(test(c).fp)*.5 || ...
        prod(max(test(c).fp,[],2)-min(test(c).fp,[],2))/prod(max(test(c).fp0,[],2)-min(test(c).fp0,[],2)) > 3 || ...
        prod(max(test(c).fp,[],2)-min(test(c).fp,[],2))/prod(max(test(c).fp0,[],2)-min(test(c).fp0,[],2)) < 1/3

            fp_err = [fp_err c];
    end
end

if verb && ~isempty(fp_err)
    fh = figure;
    set(fh,'Name','Misdetected test faces.')
    for c = fp_err
        img = imread([test(c.path),test(c).name]);
        imshow(img,[])
        hold on
        rectangle('Position',[test(c).rect(1),test(c).rect(2),test(c).rect(3),test(c).rect(4)],'EdgeColor','red')
        plot(test(c).fp(1,:),test(c).fp(2,:),'y+')
        hold off
        drawnow
        %pause
    end
end
test(fp_err) = [];
fprintf('\n Removed %d misdetected test images.\n',length(fp_err))

%% SDM on test images

% fp0_tmp = cellfun(@(x) x(1), {test.fp0});
% [test(:).fp0] = deal(fp0_tmp{:});
[test.fpt] = test.fp0;

if verb, figure; end

for c = 1 : length(test)

    img = imread([test(c).path,test(c).name]);
    if ~isa(img,'uint8')
        img = im2uint8(img);
    end
    img = uint8(normalize8(mean(img,3)));
    
    for k = 1 : niter 
           
        kpt_scale =  round(kscale(k)*mean([test(c).rect(3) test(c).rect(4)]));
        
        if k == -1 % k == 1
            % Pre-initialization
            [xi,yi] = meshgrid(linspace(test(c).rect(1),test(c).rect(1)+test(c).rect(3),10),...
                               linspace(test(c).rect(2),test(c).rect(2)+test(c).rect(4),10));
            pts = [xi(:),yi(:)];
            %kpt_scale =  round(kscale(1)*mean([test(c).rect(3) test(c).rect(4)]));
            %[fks_pi3,~] = extractHOGFeatures(img,pts,'CellSize',[kpt_scale kpt_scale],'BlockSize',[blockSize blockSize]);
            kpt_scale_pi =  round(kscale(1)*mean([test(c).rect(3) test(c).rect(4)]));
            cvkpt = cell2struct([num2cell(pts,2), repmat(num2cell(kpt_scale_pi),size(pts,1),1)], [{'pt'};{'size'}], 2); 
            fks_pi3 = extractor.compute(img, cvkpt);
            
            if size(fks_pi3,1) < length(pts)
                border_size = 100;
                img_tmp = [zeros(size(img,1),border_size),img,zeros(size(img,1),border_size)];
                img_tmp = [zeros(border_size,size(img_tmp,2));img_tmp;zeros(border_size,size(img_tmp,2))];
                %[fks_pi3,~] = extractHOGFeatures(img_tmp,pts+border_size,'CellSize',[kpt_scale kpt_scale],'BlockSize',[blockSize blockSize]);
                cvkpt = cell2struct([num2cell(pts+border_size,2), repmat(num2cell(kpt_scale_pi),size(pts,1),1)], [{'pt'};{'size'}], 2); 
                fks_pi3 = extractor.compute(img_tmp, cvkpt);
                
            end
            fks_pi3 = fks_pi3(:)' - fks_pi2_mean;
%             fks_pi3 = coeff'*fks_pi3;
            %param_pi3 = R_pi'*fks_pi3;% + param_pi_mean;
            param_pi3 = fks_pi3*R_pi + param_pi_mean;
%             imshow(img)
%             hold on
%             plot(test(c).fpt(1,:),test(c).fpt(2,:),'r.')
          
            test_pts_mean = repmat(mean(test(c).fpt,2),1,size(test(c).fpt,2));
            test(c).fpt = param_pi3(1)*(test(c).fpt - test_pts_mean) + repmat(param_pi3(2:3)'.*test(c).rect(3:4)',1,size(test(c).fpt,2)) + test_pts_mean;

%             plot(test(c).fpt(1,:),test(c).fpt(2,:),'b.')
%             plot(test(c).fp(1,:),test(c).fp(2,:),'g.')
%             drawnow; hold off; pause
      
        end
        
        %[fks,~] = extractHOGFeatures(img,test(c).fpt','CellSize',[kpt_scale kpt_scale],'BlockSize',[blockSize blockSize]);
        cvkpt = cell2struct([num2cell(test(c).fpt',2), repmat(num2cell(kpt_scale),size(test(c).fpt',1),1)], [{'pt'};{'size'}], 2); 
        fks = extractor.compute(img, cvkpt);
        if fks < size(data_avg,2)
            border_size = 150;
            img_tmp = [zeros(size(img,1),border_size),img,zeros(size(img,1),border_size)];
            img_tmp = [zeros(border_size,size(img_tmp,2));img_tmp;zeros(border_size,size(img_tmp,2))];
            %[fks,~] = extractHOGFeatures(img_tmp,test(c).fpt'+border_size,'CellSize',[kpt_scale kpt_scale],'BlockSize',[blockSize blockSize]);
            cvkpt = cell2struct([num2cell(test(c).fpt'+border_size,2), repmat(num2cell(kpt_scale),size(test(c).fpt',1),1)], [{'pt'};{'size'}], 2); 
            fks = extractor.compute(img_tmp, cvkpt);
        end
        fks = (fks(:) - fk_mean_s{k})';
        
        if 1
            tmp = test(c).fpt(:)' + fks*Rb{k} + dxk_mean{k};
        else
            tmp = test(c).fpt(:)' + (Y{k}'*((K{k}+lambda(k)*eye(size(Y{k},1)))\kernelmatrix('poly',X{k}',fks',0,3)))' + dxk_mean{k};
        end
        
        test(c).fpt = reshape(tmp,2,[]);
    end
  
    if verb
        imshow(img)
        hold on
        plot(test(c).fpt(1,:),test(c).fpt(2,:),'c.')
%         plot(test(c).fp(1,:),test(c).fp(2,:),'g.')
%         plot(test(c).fp0(1,:),test(c).fp0(2,:),'r.')
        drawnow; hold off; %pause
        
    end  

    fprintf('Test image #%d/%d localized.\n',c,length(test))
   
end

%% Cumulative error distribution curve

me = [];
for c = 1 : length(test)
    fp_ = [test(c).fp]; 
    fpt_ = [test(c).fpt];
    eye_dist = sqrt(sum((mean(fp_(:,37:42))-mean(fpt_(:,43:48))).^2));
    err_dist = sqrt(sum((fp_-fpt_).^2));
    me(c) = mean(err_dist)/eye_dist;
end

xe = linspace(min(me),max(me),100);
ye(size(xe)) = single(0);
for c = 1 : length(xe)
    ye(c) = sum(me<=xe(c))/length(me);
end

figure
plot(xe,ye,'-r')
xlabel('Distance metric')
ylabel('Fraction of Images')
title('Cumulative Error Distribution')
axis([0 .5 0 1])
