% Average face landmarks

dir = 'I:\detection\300w\annotations\lfpw\trainset\';
datad = cellstr(ls([dir,'*.pts']));

for c = 1 : length(datad)
    datas(c).fp = ptsread([dir,datad{c}])';
end

as = [];
for c = 1 : length(datas)
    as1 = [datas(c).fp(1,:) - min(datas(c).fp(1,:)); datas(c).fp(2,:)];
    as2 = [max(datas(c).fp(1,:)) - datas(c).fp(1,:); datas(c).fp(2,:)];
    as2 = [as2(:,17:-1:1) as2(:,27:-1:18) as2(:,28:31) as2(:,36:-1:32) as2(:,46:-1:43) as2(:,48:-1:47) as2(:,40:-1:37) as2(:,42:-1:41) as2(:,55:-1:49) as2(:,60:-1:56) as2(:,65:-1:61) as2(:,68:-1:66)];
    as(:,:,c) =  mean(cat(3,as1,as2),3);
end
as = round(mean(as,3));
lm_avg = [as(1,:)-min(as(1,:))+1; as(2,:)-min(as(2,:))+1]';
