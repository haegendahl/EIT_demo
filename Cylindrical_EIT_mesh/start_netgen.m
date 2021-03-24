function [ok] = start_netgen(geofile,meshfile,ng_params)

ok = 0;
ngsize = size(ng_params);

if nargin<3, error('not enough parameters'), end
if ~iscell(ng_params), error('ng_params not a cell-array!'), end
if prod(ngsize)~=max(ngsize), error('ng_params size mismatch!'), end
if mod(max(ngsize),2), error('ng_params length not even!'), end


%% write the options file for Netgen %%
ngopt_fid = fopen('ng.opt','w');
if ngopt_fid==-1, error('couldn''t open options file: ng.opt!'), end

for ii=1:2:max(ngsize)-1
  ngopt_str = parse_ngoptions(ng_params{ii},ng_params{ii+1});
  fprintf(ngopt_fid,'%s\n',ngopt_str);
end
fclose(ngopt_fid);


%% start Netgen %%
finelevel = '';

% Windows 
%ngrun_str = sprintf('"c:\\program Files (x86)\\netgen-4.9.13_Win32\\bin\\netgen.exe" %s -batchmode -geofile=%s  -meshfile=%s ',finelevel,geofile,meshfile);

%linux command line
ngrun_str = sprintf('netgen %s -batchmode -geofile=%s  -meshfile=%s ', ...
         finelevel,geofile,meshfile);
disp(ngrun_str);

ng_status = system(ngrun_str);

%if ng_status~=0; error('something went wrong?!'), end

ok = 1;
return



%% Parse function %%
function [str] = parse_ngoptions(opt,val)

str = '';

switch upper(opt)
 case 'MESH GRANULARITY'
  if ~isnumeric(val), error('Mesh granularity value is not numeric!'), end
  val = uint8(fix(val));
  str = ['meshoptions.fineness ' num2str(val)];
  
 case '2ND ORDER ELEMENTS'
  if (~isnumeric(val)) | (~(val==0 | val==1)), 
    error('2nd order element value should be 0 or 1!')
  end
  val = uint8(fix(val));
  str = ['options.secondorder ' int2str(val)];
  
 case 'MAX MESH-SIZE'
  if ~isnumeric(val), error('Max mesh-size value is not numeric!'), end
  str = ['options.meshsize ' num2str(val)];
  
 case 'MESH-SIZE GRADING'
  if ~isnumeric(val), error('Mesh-size grading value is not numeric!'), end
  str = ['options.grading ' num2str(val)];
  
 case 'ELEMENTS PER CURVATURE'
  if ~isnumeric(val), error('Elements per curvature value is not numeric!'), end  
  str = ['options.curvaturesafety ' num2str(val)];
  
 case 'ELEMENTS PER EDGE'
  if ~isnumeric(val), error('Elements per edge value is not numeric!'), end  
  str = ['options.segmentsperedge ' num2str(val)];
  
 otherwise
  error(sprintf('Unknown meshing parameter: %s',opt))
end
