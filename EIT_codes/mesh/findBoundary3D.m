function [bndnodes,bndfaces] = findBoundary3D(H,g);
% finds the nodes (and faces) of 3D triangular mesh.
% Uses surftri.m by Per-Olof Persson

bndfaces = surftri(g,H);
bndnodes = unique(bndfaces(:));
