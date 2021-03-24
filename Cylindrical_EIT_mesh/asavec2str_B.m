function [boxstr] = asavec2str_B(g,Nel,fname)
% electrode placement for netgen. output is a text string containing netgen code for
% the electrodes (polyhedrons).
%   input:     g = corner point coordinates for polyhedrons (electrodes)
%            Nel = # of electrodes
%          fname = filename (optional)

boxstr = cell(Nel,1);
for nel=1:Nel
  vecid = 1 + 8*(nel-1);
  pikkug=g(vecid:vecid+7,:);
  str = ['solid box' num2str(nel) ' = polyhedron('];
  for ii=1:8,
     for jj=1:2,
      
      str = [str num2str(pikkug(ii,jj)) ','];
     end
     str = [str num2str(pikkug(ii,3))];
     str= [str ';'];
  end
  str = [str '; 1, 5, 8; 1, 8, 4; 4, 8, 7; 4, 7, 3; 3, 6, 2; 3, 7, 6; 2, 5, 1; 2, 6, 5; 1, 4, 3; 1, 3, 2; 6, 7, 5; 7, 8, 5);'];
  boxstr{nel} = str;
end

% write to file
if nargin==3
  fid = fopen(fname,'wt');  % open for writing (text mode)
  if fid==-1
    error('Cannot open file: %s',fname)
  else
    for ii=1:Nel
      fprintf(fid,'%s\n',boxstr{ii});
    end
    fclose(fid);
  end
end
