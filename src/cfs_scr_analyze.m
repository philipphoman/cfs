function out = analyze_cfs_new(in)
% Analyze CFS data
%
% Algorithm:
% For each subject:
% 1. Load .acq raw file with acq_load
% 2. Save file as .mat file
% 3. Import data with pspm, maker sure import.flank ='ascending'
% 4. Trim data, from 3. marker + 20 sec until 31. marker + 20 sec
% 5. Set up timing file: events matrix, use data.data 
% 6. Set up first level: dcm, 4 conditions (CS+ CS- CS+- CS-+)
% 7. Define contrasts: CS+ > CS-, CS-+ > CS+-
% End For
% 8. Define 2nd level: load each dcm, one-sample ttest
% 9. Export statistics
%
% NB: dcm.stats will contain the trialwise estimates. 
% The columns correspond to dcm.names, i.e., usually the amplitudes
% are in the first column.
%
% PH, 20151203

% root directory
out = 1;
rootdr = pwd;

% do preprocessing?
dopreproc = 1;

% condition length
epochlength = 6;

% marker stuff for trimming
from = 14;
to = 21;
refmarker1 = [3, 35];
refmarker2 = [4, 36];

f=dir(fullfile(rootdr,'subjects/*acq'));

if nargin==1
  str=in;
else
  e=numel(f);
  str=['1:' num2str(e)];
end
nn=regexp(str,':','split');
n=[str2num(nn{1}):str2num(nn{2})];
for i=1:numel(n)
  %[a1,a2,a3] = fileparts(f(n(i)).name);
  %fp = regexp(a2,'-','split');
  %str2=[str2 ' -i natives/' f(i).name];
  %if isdir(fullfile(dr,'subjects',fp{1}))
  %if exist(fullfile(dr,'subjects',fp{1},'mri','orig','001.mgz'))
  %end
  subjfiles{i}=fullfile(rootdr,'subjects',f(n(i)).name);
end

dbadvice = 1;
dooverwrite = 0;

%for grp = 1:numel(group)
  for subj = 1:numel(subjfiles)
    
    if dopreproc == 1

      % check if dcm file is present
      s = subjfiles{subj};
      [p,f,x] = fileparts(s);
      df = fullfile(p,['tscr_' f '_dcm.mat']);
      if exist(df)
        fprintf('Found dcm file %s, exiting...\n',df);
        continue
        clear df p f x s
      end
    
      % make .mat-file from .acq-file
      s = subjfiles{subj};
      m = load_acq(s);
      data = m.data;
      [p,f,x] = fileparts(s);
      save(fullfile(p,f),'data');
      subjfiles{subj} = fullfile(p,[f '.mat']);
      s = subjfiles{subj};
      clear p f x
      
      % import file with pspm
      import{1} = struct('type','scr','channel',1,'sr',200,'transfer','none');
      import{2} = struct('type','marker','channel',2,'sr',200,'flank','ascending');
      %options.overwrite = dooverwrite;
      options.overwrite = 1;
      scr_import(s,'mat',import,options);
      [p,f,x] = fileparts(s);
      subjfiles{subj} = fullfile(p,['scr_' f '.mat']);
      s = subjfiles{subj};
      clear p f x
      
      % trim file with pspm
      % determine whether 35 (quasi-default) or 36 markers are present
      d = load(s);
      if size(d.data{2}.markerinfo.value,1) == 35
        refmarker = refmarker1;
        scr_trim(s,from,to,refmarker,struct('overwrite',1,'marker_chan_num',0));
      elseif size(d.data{2}.markerinfo.value,1) == 36
        refmarker = refmarker2;
        scr_trim(s,from,to,refmarker,struct('overwrite',1,'marker_chan_num',0));
      else
        % cut away the first four markers and then analyze the rest 
        refmarker = 'marker';
        scr_trim(s,60,15,refmarker,struct('overwrite',1,'marker_chan_num',0));
      end
      %scr_trim(s,from,to,refmarker,struct('overwrite',1,'marker_chan_num',0));
      [p,f,x] = fileparts(s);
      subjfiles{subj} = fullfile(p,['t' f '.mat']);
      s = subjfiles{subj};
      
      % set up timing file
      load(s);

      % according to Dominik, it might be good to model timings differently given 
      % the long SOA use here (see his email from 8/4/16
      dc = 0;
      if dbadvice == 1
        for dd=1:size(data{2}.data,1)
          %dc = dc + 1;
          %dtmp = data{2}.data(dd);
          %aevents(dd,1) = dtmp;
          %eevents(dc,1) = dtmp;
          %dc = dc + 1;
          %aevents(dc,1) = dtmp;
          %eevents(dc,1) = dtmp + 4;
          %dc = dc + 1;
          %aevents(dc,1) = dtmp;
          %eevents(dc,1) = dtmp + 6;
        end
        %events{1} = [data{2}.data data{2}.data+epochlength];
        events{1} = [data{2}.data];
        events{2} = [data{2}.data + 4];
        events{3} = [data{2}.data + epochlength];
      else
        events{1} = [data{2}.data data{2}.data+epochlength];
        events{2} = [data{2}.data+epochlength];
      end
      clear dc
      %events{3} = [data{2}.data(17:32) data{2}.data(17:32)+epochlength];
      %events{4} = [data{2}.data(17:32)+epochlength];
      timingfn = fullfile(p,['t' f '_timings.mat']);
      save(timingfn,'events');
      clear p f x events
    end      

    % set up dcm first level with pspm
    % reload s again
    s = subjfiles{subj};
    [p,f,x] = fileparts(s);
    %subjfiles{subj} = fullfile(p,['tscr_' f '.mat']);
    subjfiles{subj} = fullfile(p,[f '.mat']);
    s = subjfiles{subj};

    % reload timing file
    %timingfn = fullfile(p,['tscr_' f '_timings.mat']);
    timingfn = fullfile(p,[f '_timings.mat']);

    %model.modelfile = fullfile(p,['tscr_' f '_dcm.mat']);
    model.modelfile = fullfile(p,[f '_dcm.mat']);
    model.datafile = s;
    model.timing = timingfn;
    model.norm = 1; % normalize model (recommended in within-subject designs)
    options = struct(...
      'crfupdate',0,...
      'indrf',0,...
      'getrf',0,...
      'rf',0,...
      'depth',2,...
      'sfpre',2,...
      'sfpost',5,...
      'sffreq',0.5,...
      'sclpre',2,...
      'sclpost',5,...
      'aSCR_sigma_offset',0.3,...
      'dispwin',0,...
      'overwrite',dooverwrite,...
      'dispsmallwin',0 ...
    );

    dcm = scr_dcm(model,options);
    %dcmfiles{grp}{subj} = model.modelfile;
    
    % set up contrasts with pspm
    %scr_con1(model.modelfile, {'CS++ > CS--','CS-+ > CS+-'},...
    %  eval(['convec.' exporder{grp}{subj}]),'param',0,struct('zscored',0));   
    
  end
  
  % group level
  clear options
  %outfile = fullfile(rootdr,['dcm_grouplevel_' num2str(grp)]);
  %con = 'all';
  %options.overwrite = 0;
  %connames = 'number';
  %scr_con2(dcmfiles{grp},outfile,con,connames,options);
    
end  
