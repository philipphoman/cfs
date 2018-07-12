rootdr=pwd;
subjdr=fullfile(rootdr,'subjects');
fn=dir(fullfile(subjdr,'tscr*dcm.mat'));
M=[];
for i=1:numel(fn)
  [p f x]=fileparts(fn(i).name);
  fname=fullfile(subjdr,[f x]);
  subid=f(6:9);
  subidnum=str2num(subid);
  fprintf('Processing subject %i (%s)...',subidnum,subid);
  d=load(fname);
  dcm{i}=d.dcm.stats;
  try
    M=[M; repmat(subidnum,size(d.dcm.stats,1),1) [1:size(d.dcm.stats,1)]' d.dcm.stats];
    fprintf('done\n');
  catch ME
    fprintf('failed with error %s\n',ME.message);
  end
end