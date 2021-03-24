function slices = GetSlices(g,theta)


slice_index = unique(g(:,3));

slice_index = slice_index(2:end-1);

for ii = 1:numel(slice_index)
    
    ind = g(:,3)== slice_index(ii);
    slices(ii).g = [g(ind,1), g(ind,2)];
    slices(ii).H = delaunay(g(ind,1),g(ind,2)); 
    slices(ii).theta = theta(ind,end);
    
end

