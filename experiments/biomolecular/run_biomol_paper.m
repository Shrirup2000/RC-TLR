function out = run_biomol_paper(dataRoot,options)
%RUN_BIOMOL_PAPER Run the standalone RC-TLR experimental suite.
% out=run_biomol_paper(); % uses the module's data directory
% out=run_biomol_paper(dataRoot,struct('qrTol',4e-8,'target',1e-7));
    here=fileparts(mfilename('fullpath')); addpath(here,'-begin');
    if nargin<1 || isempty(dataRoot), dataRoot=fullfile(here,'data'); end
    if nargin<2 || isempty(options), options=struct; end
    assert(isfolder(dataRoot),'Data folder does not exist: %s',dataRoot);
    plan=biomol_experiment_plan(options); names=unique({plan.filename});
    fprintf('RC-TLR experiment suite: %d configurations.\n',numel(plan));
    restore_biomol_data(dataRoot,names);
    files=dir(fullfile(dataRoot,'**','*.mat'));
    problems=strings(0,1);
    for k=1:numel(names)
        count=nnz(strcmp({files.name},names{k}));
        if count~=1
            problems(end+1)=sprintf('%s: found %d copies; exactly one is required.',names{k},count); %#ok<AGROW>
        end
    end
    assert(isempty(problems),'Input check failed:\n%s',strjoin(problems,newline));
    out=rc_tlr_biomol_experiment(dataRoot,options);
    if ~out.AllRequestedCasesOK || strlength(out.PlotError)>0
        warning('biomol:Incomplete','Some requested cases or plots are incomplete. Inspect all_cases.csv and *_FAILED.txt in %s.',out.Settings.outDir);
    end
end
