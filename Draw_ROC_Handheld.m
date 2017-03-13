

X1 = load('Handheld_Front_Vito_X.mat');
Y1 = load('Handheld_Front_Vito_Y.mat');
X2 = load('Handheld_Front_2D_X.mat');
Y2 = load('Handheld_Front_2D_Y.mat');
X3 = load('Handheld_Front_PreVGG_X.mat');
Y3 = load('Handheld_Front_PreVGG_Y.mat');
X4 = load('Handheld_Front_Orig_X.mat');
Y4 = load('Handheld_Front_Orig_Y.mat');

w = 1.4;

plot(X1.X1,Y1.Y1,'g','LineWidth',w);
hold on
plot(X2.X2,Y2.Y2,'--k','LineWidth',w);
hold on
plot(X3.X3,Y3.Y3,'-.k','LineWidth',w);
hold on
plot(X4.X4,Y4.Y4,':k','LineWidth',w);

xlabel('False Accept Rate');
ylabel('True Accept Rate');
%xticklabels({'top-1','top-5','top-10'});
h = legend('Our method (asymmetric)','2D aligned','Pre-trained (2D aligned)','No pre-processing','Location','southeast');
set(h,'FontSize',7.5);  
%set(gca, 'XTickLabel', {'top-1','top-5','top-10'});
set(gcf, 'papersize', [5 5]);
   set(gcf, 'paperposition', [0 0 5 5]);
   print -dpdf C:\Users\Sandipan\Desktop\NDFaceNet_New\FG2017_latex_template\FG2017_latex_template\Images\ROC_Handheld_Front_OldVito_asym.pdf

% val = 0.02; %value to find
% tmp = abs(X1-val);
% [idx idx] = min(tmp); %index of closest value
% %closest = X1(idx) %closest value
% Y1(idx)
% tmp = abs(X2-val);
% [idx idx] = min(tmp); %index of closest value
% %closest = X2(idx) %closest value
% Y2(idx)
% tmp = abs(X3-val);
% [idx idx] = min(tmp); %index of closest value
% Y3(idx)
% tmp = abs(X4-val);
% [idx idx] = min(tmp); %index of closest value
% Y4(idx)

