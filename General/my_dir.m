function files = my_dir(path_name, suffix, begin)
% files = MY_DIR(path_name,begin,suffix) returns files and directories of   
% directory $path_name. They will be sorted by name.
% If begin is given MY_DIR returns only files/directories beginnig with 
% $begin.
% If $suffix is given MY_DIR returns only files with ending of $suffix
% $suffix might be '*'.

if(nargin < 1),
    path_name ='.';
end,
if(nargin < 2),
    tmp = dir(path_name);
    start = 3;
else
    tmp = dir([path_name, '\*.', suffix]);
    start = 1;
end,
if(nargin < 3),
    begin = [];
end,
str = ''; m=0;

% get the [begin *]-files only
for (n=start:length(tmp)),
   tmp2 = getfield(tmp(n),'name');
   if length(tmp2)>=length(begin),
      if(isempty(begin) | strcmp(tmp2(1:length(begin)),begin)) , 
         m=m+1;
         tmp3{m} = tmp2;
      end,
   end,
end,
% sort file names (only possible as matrix)(tmp3 is cell array!)
if (exist('tmp3')),
    for (n=1:length(tmp3)),
         files(n,1:length(tmp3{n})) = tmp3{n};
    end,
    files = sortrows(files);
else,
    files = [];
end,    
    