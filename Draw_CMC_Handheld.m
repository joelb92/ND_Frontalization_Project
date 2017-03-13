y1 = [54.24,60.17,64.28,66.81,68.91,70.48,71.35,72.58,73.62,74.85,75.28,75.98,76.59,77.47,77.73]; %Vito
y2 = [13.54,17.21,20.52,22.27,24.02,25.85,26.81,28.47,29.52,30.83,31.97,32.66,33.71,34.67,35.46]; %2D
y3 = [28.82,35.11,37.47,40.0,41.92,43.93,45.59,46.81,47.77,48.82,49.96,50.92,52.05,52.66,53.36]; %Pre
y4 = [3.67,4.98,6.29,7.6,9.34,10.22,11.0,11.97,12.84,13.28,14.32,15.11,15.63,16.42,16.68]; %Orig
x = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
w = 1.4;
plot(x,y1,'g','LineWidth',w);
hold on
plot(x,y2,'--k','LineWidth',w);
hold on
plot(x,y3,'-.k','LineWidth',w);
hold on
plot(x,y4,':k','LineWidth',w);

xlabel('Rank');
ylabel('Recognition Accuracy (%)');
%xticklabels({'top-1','top-5','top-10'});
h = legend('Our method (asymmetric)','2D aligned','Pre-trained (2D aligned)','No pre-processing','Location','southeast');
set(h,'FontSize',7.5);  
%set(gca, 'XTickLabel', {'top-1','top-5','top-10'});
set(gcf, 'papersize', [5 5]);
   set(gcf, 'paperposition', [0 0 5 5]);
   print -dpdf C:\Users\Sandipan\Desktop\NDFaceNet_New\FG2017_latex_template\FG2017_latex_template\Images\CMC_Handheld_OldVito_Front_Newest.pdf
