%clear all;
%clc;
trueMatchScores2 = [];
falseMatchScores2 = [];
FeatFile = importdata('C:\Users\Sandipan\Desktop\NDFaceNet_New\PaSC_Old_Vito\Full_handheld_Feat_Orig_Averaged.txt');
FeatFile1 = importdata('C:\Users\Sandipan\Desktop\NDFaceNet_New\PaSC_Old_Vito\Full_handheld_Feat_Orig_Averaged.txt');
for i=1:length(FeatFile.textdata)
    i
    currentFile1 = FeatFile.textdata{i};
    currentLabel1 = currentFile1(1:6);
    currentFeat1 = FeatFile.data(i,:);
    for j = 1:length(FeatFile1.textdata)
        currentFile2 = FeatFile1.textdata{j};
        currentLabel2 = currentFile2(1:6);
        currentFeat2 = FeatFile1.data(j,:);
        %score = norm(currentFeat1 - currentFeat2);
        score = pdist2(currentFeat1,currentFeat2,'cosine');
        if strcmp(currentLabel1, currentLabel2)
            trueMatchScores2(end+1) = score;
        else
            falseMatchScores2(end+1) = score;
        end
    end
end

scores2 = [trueMatchScores2 falseMatchScores2; ones(1,length(trueMatchScores2)) zeros(1,length(falseMatchScores2))];
[X4,Y4] = perfcurve(scores2(2,:), scores2(1,:),0);

