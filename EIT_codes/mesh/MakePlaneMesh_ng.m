function  MakePlaneMesh_ng(dims,elcenters,eldims,eltype,meshname)

xlen = dims(1);
ylen = dims(2);
zlen = dims(3);

Nel = size(elcenters,1);

%% here goes...
writestr = ['algebraic3d;' char([10 10]) 'solid top = plane(0,0,'];
writestr = [writestr num2str(zlen) ';0,0,1);' char(10)];
writestr = [writestr 'solid box = plane(0,0,0;0,0,-1)' char([10 9]) 'and plane(0,0,' ...
                    '0;-1,0,0)' char([10 9]) 'and plane(0,0,0;0,-1,0)' char([10 9])];
writestr = [writestr 'and plane(' num2str(xlen) ',0,0;1,0,0)' char([10 9])];
writestr = [writestr 'and plane(0,' num2str(ylen) ',0;0,1,0)' char([10 ...
                    9]) 'and top;' char([10 10])];

if strcmp(eltype,'rect')
  disp('Not implemented yet!')
  return

      
elseif strcmp(eltype,'circ')
  
  writestr2 = '';
  for ii=1:Nel
    cx = elcenters(ii,1);
    cy = elcenters(ii,2);
    radius = eldims(ii,1);
    xystr = [num2str(cx) ',' num2str(cy) ','];
    writestr = [writestr 'solid el' num2str(ii) ' = top' char([10 9])];
    writestr = [writestr 'and plane(0,0,0;0,0,-1)' char([10 9])];
    writestr = [writestr 'and cylinder(' xystr '0;' xystr num2str(zlen) ...
                ';' num2str(radius) ');' char(10)];    
    
    writestr2 = [writestr2 'tlo el' num2str(ii) ' top -col=[1,0,0];' char(10)];
  end
  writestr = [writestr char(10) 'tlo box -transparent;' char(10) writestr2];

elseif strcmp(eltype,'inverse')
  H2d = delaunay(g2d(:,1),g2d(:,2));
  elindall = [];
  
else
  disp('error!')  
  g2d = 0;
  return
end

fid = fopen(meshname,'w');
if fid~=-1
  fprintf(fid,'%s\n',writestr);
  fclose(fid);
else
  error(['couldn''t open file: ' meshname])
end

  
