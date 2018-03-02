function [imgOut,roi] = maskFaceByPoints(I_Q,fidu_XY)
        xcords = double(fidu_XY(:,1));
        ycords = double(fidu_XY(:,2));
        hull = convhull(xcords,ycords);
        xcords = xcords(hull);
        ycords = ycords(hull);
        mask = repmat(uint8(poly2mask(xcords,ycords,size(I_Q,1),size(I_Q,2))),[1,1,3]);
        imgOut = (I_Q).*(mask);
        roi = [min(xcords) min(ycords) max(xcords)-min(xcords) max(ycords)-min(ycords)];
end