function [bndnodes,bndedges] = findBoundary2D(H,g)
% finds boundary nodes (and edges) of 2d triangular mesh


H = uint32(H);

% get edges, non-dublicates are surface edges
edges = sort([H(:,[1 2]);H(:,[2 3]);H(:,[3 1])],2);

[useless I J] = unique(edges,'rows');
count = histc(J,1:max(J));
nonid = find(count==1);

bndedges = edges(I(nonid),:);
bndnodes = unique(bndedges(:));
