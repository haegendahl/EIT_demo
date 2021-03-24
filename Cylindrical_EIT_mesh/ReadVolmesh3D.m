function [H,g,h,bnd] = ReadVolmesh3D(volfile)
% READVOLMESH3D Read Netgen's 3D vol-files.
%    [TET,NODE,TRI,BND] = READVOLMESH3D Reads a 3D mesh file in
%    Netgen's VOL-format. A file chooser dialog is opened for
%    choosing the file. The output parameters are as follows:
%       TET .... matrix index array for volume elements (tetrahedrons)
%       NODE ... n-by-3 matrix array of the nodal points (x,y,z)
%       TRI .... matrix index array for surface triangulation
%       BND .... m-by-1 array of boundary indices for TRI
%
%    [TET,NODE,TRI,BND] = READVOLMESH3D(FILE) Works as above but
%    the volmesh-file can be specified directly in char array
%    FILE instead of opening a file dialog.
%
%    This version supports only tetrahedral meshes and 1st or 2nd
%    order basis for TET and TRI.
%
%
% Copyright University of Eastern Finland, 
% Department of Applied Physics
% P.O. Box 1627, FIN-70211 Kuopio, Finland

% K. Karhunen 2011


global DEBUG;
DEBUG = false;


% open a file dialog (if no filename given)
if nargin==0
  [fname,fpath] = uigetfile('*.vol','Select a vol-file');
  if ~ischar(fname)
    error('Cancel!')
  end
  
  volfile = fullfile(fpath,fname);
end

% open the file for reading (if possible)
fid = fopen(volfile,'r');
if fid==-1
  error(sprintf('Couldn''t open file: %s',volfile))
end


% read header
fseek(fid,0,-1);   % jump to beginning of file (forced)

sline = fgetl(fid);
if ~strcmp(sline,'mesh3d')
  myerror(fid,'Unknown mesh format!')
end

sline = fgetl(fid);
if strcmp(sline,'dimension')
  sline = fgetl(fid);
  if ~strcmp(sline,'3')
    myerror(fid,'Mesh dimension not 3!')
  end
else
  myerror(fid,'Mesh dimension tag not found!')
end

%% header checks passed... start reading the data %%



% define subheaders to be read
subhead = {'surfaceelements',...    % surface elements, Netgen version >=4.9.x
           'surfaceelementsgi',...  % surface elements, version <4.9.x
           'volumeelements',...
           'points',...
           'edgesegmentsgi2'};
          

ns = 4;  % # of subheaders needed to be found

% automatic space generation...
spc = aaddspace(subhead);

% loop until all read
cnt = 0;
while 1
  sline = fgetl(fid);
  
   % check for an error, end of file or all done
  if (~ischar(sline)) || (cnt==ns), break, end 
  
  test = 0;
  for ii=1:length(subhead)
    if strcmp(sline,subhead{ii})
      test = ii;
      cnt  = cnt + 1;
      fprintf(' %s%s: ',sline,spc{ii})
      break;
    end
  end
  
  if test
    nel = str2num(fgetl(fid));  % should be # of elements
    if isempty(nel)
      myerror(fid,'file corrupt!')
    end
    
    % try reading a whole data block at once
    data{cnt,1} = readvoldata(fid,nel,test);
    data{cnt,2} = test;
    fprintf('%d\n',max(size(data{cnt,1})))
  end
  % nop
end

fid = fclose(fid);
if fid<0
  disp('Warning: couldn''t close the vol-file!')
end


% parse data for output (edge segments neglected)
[H,h,g,bnd] = parsevoldata(data);
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% sub function to parse raw voldata
function [H,h,g,b,e] = parsevoldata(data)

nd = size(data,1);  % # of data fields

% initialize output
H = [];
h = [];
g = [];
b = [];
e = cell(1,1);

% loop through
foo = 1;
for ii=1:nd
  type = data{ii,2};
  dlen = size(data{ii,1},2);
  
  if (type==1) || (type==2)
    % surface mesh (faces)
    h = data{ii,1}(:,6:dlen);  % only face indices
    b = data{ii,1}(:,2);       % boundary data
  elseif type==3
    % volume mesh
    H = data{ii,1}(:,3:dlen);
  elseif type==4
    % node coordinates
    g = data{ii,1};
  else
    % edge segments + other stuff if it exists   
    e{foo} = data{ii,1};
    foo    = foo + 1;
  end
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% sub function for reading data blocks
function [data] = readvoldata(fid,nlines,id);
global DEBUG

% read two lines to check data consistency (include LF)
sline1 = fgets(fid);  
sline2 = fgets(fid);

% determine number of bytes to be read
% this may need further updates as Netgen versions change
if length(sline1) == length(sline2)
  % line terminates to one LF
  nlines = nlines - 2;      % # of lines/elements
  lbytes = length(sline1);  % % of bytes on a line
elseif (max(size(sline2))==1) && (uint8(sline2(1))==10)
  % line terminates to two LF's  
  nlines = nlines - 1;
  lbytes = length(sline1) + 1;
  sline1 = [sline1 ' '];
  sline2 = [];
elseif id==5
  % edge segments (old format (messy))
  nlines = nlines - 2;
  lbytes = length(sline1) + length(sline2);
  sline1 = [sline1 sline2];
  sline2 = [];
  %keyboard
else  
  error('Data block/file corrupt!')
end

% read from file
bytes   = nlines*lbytes;
bindata = fread(fid,bytes,'uchar');


ndata = max(size(bindata));
if ndata ~= bytes
  disp('Warning: not all data bytes read!')
  disp('  (maybe partial data block/file corruption (or a bUg o_O)!)')
end

nlines = fix(ndata/lbytes);  % select only complete data lines

% reshape to numerical array
data = [sline1;sline2;reshape(bindata,lbytes,nlines)'];
data = str2num(char(data));

if DEBUG
  disp(sprintf('%d bytes read == %d rows - %d columns!',ndata+length(sline1)+length(sline2), ...
               size(data,1),size(data,2)))
end
return


% automatic add space
function [spc] = aaddspace(strcell)

smax = 0;
for ii=1:max(size(strcell))
  smax = max([smax, length(strcell{ii})]);
end
for ii=1:max(size(strcell))
  spc{ii} = char(32*ones(1,smax-length(strcell{ii})));
end
return


% error handling in a very simple (and wrong) way
% just a poor man's "check" that the file is closed, nothing else
function [foo] = myerror(fid,msg);

fclose(fid);
disp(['**  PROGRAM ERROR: ' msg ' **'])
error(' ')

% never returns
return
