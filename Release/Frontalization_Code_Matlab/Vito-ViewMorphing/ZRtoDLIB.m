function reordered = ZRtoDLIB(points)
reordered = zeros([68,2]);
% Jaw
reordered(1:9,:) = flipud(points(52:60,:));
reordered(10:17,:) = points(61:68,:);
% left eyebrow
reordered(18:22,:) = points(16:20,:);
% right eyebrow
reordered(23:27,:) = points(27:31,:);
% nose
reordered(28:32,:) = flipud(points(5:9,:));
reordered(32:33,:) = flipud(points(2:3,:));
reordered(34,:) = points(1,:);
reordered(35:36,:) = points(4:5,:);
% left eye
reordered(37:39,:) = flipud(points(13:15,:));
reordered(40:42,:) = points(10:12,:);
% right eye
reordered(43:45,:) = points(24:26,:);
reordered(46:48,:) = flipud(points(21:23,:));
% outer mouth 
reordered(49:52,:) = flipud(points(32:35,:));
reordered(53:55,:) = flipud(points(39:41,:));
reordered(56,:) = points(44,:);
reordered(57,:) = points(45,:);
reordered(58,:) = points(51,:);
reordered(59,:) = points(48,:);
reordered(60,:) = points(50,:);
% inner mouth
reordered(61:63,:) = points(36:38,:);
reordered(64,:) = points(42,:);
reordered(65,:) = points(43,:);
reordered(66,:) = points(45,:);
reordered(67,:) = points(47,:);
reordered(68,:) = points(49,:);
end