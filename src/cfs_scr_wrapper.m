% This is the wrapper to call the cfs-analyzer
function out = cfs_wrapper(input)
  %#function g_SCR f_SCR f_SF
  out=-1;
  if nargin==1
    analyze_cfs(input);
  else
    analyze_cfs;
  end
  out=0;
end
