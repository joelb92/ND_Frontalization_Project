function [new_im, grn, gcn] = transformImage(I, gr, gc, TM)

if isa(TM,'projective2d') || isa(TM,'affine2d')
    [x,y] = transformPointsForward(TM, gc, gr);
    [new_im,RJ] = imwarp(I,TM, 'OutputView',...
        imref2d(size(I),...
                [mean(x(1:end-4))-size(I,1)/2, mean(x(1:end-4))+size(I,1)/2],...
                [mean(y(1:end-4))-size(I,2)/2, mean(y(1:end-4))+size(I,2)/2]));
    xdata = RJ.XWorldLimits;
    ydata = RJ.YWorldLimits;
else
    [new_im, xdata, ydata] = imtransform(I, TM, 'size',round(1.5*[size(I,1) size(I,2)]));
end
w = xdata(2)-xdata(1);%+1;
h = ydata(2)-ydata(1);%+1;
scalex = size(new_im,2)/w;
scaley = size(new_im,1)/h;

% cartesian coordinates
go = [gc(:), gr(:)];

%% cartesian coordinates
if isa(TM,'projective2d') || isa(TM,'affine2d')
    NPt = transformPointsForward(TM, go);
else
    NPt = tformfwd(TM, go);
end
% convert to image coordinates
NPtImage(:,1) = NPt(:,1) - xdata(1);% + 1;
NPtImage(:,2) = NPt(:,2) - ydata(1);% + 1;
NPtImage(:,1) = NPtImage(:,1)*scalex;
NPtImage(:,2) = NPtImage(:,2)*scaley;
NPtImage = round(NPtImage);

grn = NPtImage(:,2);
gcn = NPtImage(:,1);
grn = reshape(grn, size(gr,1), size(gr,2));
gcn = reshape(gcn, size(gc,1), size(gc,2));

showimgs = false;
if showimgs
  I = uint8(I);
  new_im = uint8(new_im);

  figure, imshow(I), hold on;
  for i=1:size(go,1)
      plot(go(i,1), go(i,2), 'r*');
  end

  figure, imshow(new_im, 'XData', xdata, 'YData', ydata);
  hold on
  axis on
  for i=1:size(go,1)
      plot(NPt(i,1), NPt(i,2), 'r*')
  end

  tmpim = new_im;
  for i=1:size(NPtImage,1)
      tmpim = makept(tmpim, NPtImage(i,2), NPtImage(i,1), 'b');
  end
  figure, imshow(tmpim);
end

end


function I = makept(I, r, c, color)
mr = max(r-2,1);
mc = max(c-2,1);
maxr = min(r+2,size(I,1));
maxc = min(c+2,size(I,2));

I(mr:maxr, mc:maxc,:) = 0;

if strcmp(color, 'r')
  I(mr:maxr, mc:maxc,1) = 255;
else
  I(mr:maxr, mc:maxc,3) = 255;
end
end