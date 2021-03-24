load data






% ind = ginv(:,3)<14.5 & ginv(:,3) > 14.0;
% 
% 
% ginv_2D = [ginv(ind,1), ginv(ind,2)];
% 
% Hinv_2D = delaunay(ginv_2D(:,1),ginv_2D(:,2)); 
% 
% theta_2D = theta(ind,end);
% 
% close all
% figure(1)
% 
% patch('faces',Hinv_2D, 'vertices',ginv_2D, 'facevertexcdata',theta_2D,'facecolor','interp','edgecolor','none');
% %caxis([min(theta_2D),max(tehta_2D)])




slices = GetSlices(ginv,theta);


for i= 1 : numel(slices)
    figure(1)
    slice_num = i;

    PlotSolScaled(slices(slice_num ).g,slices(slice_num ).H,slices(slice_num ).theta,0.5,2);
    pause
end

for i= 1 : numel(slices)
    slice_num = i;
    [R Ai] = constructTVIntMat(slices(slice_num ).g,slices(slice_num ).H);
    PlotGradient(slices(slice_num ).theta, R,slices(slice_num ).g, slices(slice_num ).H, 2)
    pause
end