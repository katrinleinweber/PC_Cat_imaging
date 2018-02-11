function out = cfg_run_gzip_files(job)
% Run gzip on a list of files.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_gzip_files.m 6840 2016-07-25 12:21:25Z guillaume $

rev = '$Rev: 6840 $'; %#ok

if isempty(job.outdir) || isempty(job.outdir{1})
    for i=1:numel(job.files)
        out = gzip(job.files{i});
    end
else
    out = gzip(job.files, job.outdir{1});
end
if ~job.keep
    delete(job.files{:});
end