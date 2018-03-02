function [features, noValidPoints] = MyextractHOGFeatures1(I, points, csize, bsize, overlap, nbins, signedhist)
% The function computes HOG features from specified points. It was writen
% specifically with speed in mind. There is no paraemtere checking i hte
% mex file that is called, the only parameter checking that is done is here
% in this file. If you do not call the function as I intended, you will
% produce a segmentaion fault - I promise.
% 
% The function needs to take all seven input parameters. If youwant to
% provide less and resort to defaults, change this file and provide
% defaults.
% 
% Prototype:
%   [features, noValidPoints] = MyextractHOGFeatures(I, points, csize, bsize, overlap, nbins, signedhist)
% 
% Inputs:
% 
% Outputs:
%   features ... a NxM matrix, where N is the number of provided points
%                from which to extract the HOGs and M is the dimensionality
%                of the HOGs, which equals bsize*bsize*nbins
%   noValidPoints ... number of valid points from which the HOGs can be
%                   computed
% 
% Inputs:
%   points      ... a Nx2 matrix with coordinates of x,y coordinats of N 
%                   points 
%   csize       ... cell size in pixels, usually 8 or 16
%   bsize       ... block size in terms of how many cells there are in one 
%                   block, e.g., bsize = 2 if you want 2x2 cells in one 
%                   block 
%   overlap     ... this has to be zero 
%   nbins       ... number of bins in the HOG histograms, usually 9
%   signedhist  ... flag to use either the total range orientation range 
%                   0-360 (for true) or half the range 0-180 (for false)
% 
% 
% Author: Vitomir Struc
% Date: 31.8.2015
% Copyright: Faculty of Electrical Engineering, University of Ljubljana,
% University of Notre Dame, 2015

%% Parameter check
I = uint8(I);
if overlap~=0
    overlap=0;
end

if mod(bsize,2)~=0
    bsize=bsize+1;
end

if mod(csize,2)~=0
    csize=csize+1;
end

%% Call MEX function 
[features, noValidPoints] = vitos_hog2(I,points, csize, bsize, overlap, nbins, signedhist);





















