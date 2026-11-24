function [jobs,c] = biomol_experiment_plan(options)
%BIOMOL_EXPERIMENT_PLAN Enumerate RC-TLR studies without loading data.
% [jobs,cfg]=biomol_experiment_plan(); % 55 default configurations
% p counts nodes/coefficient slots, as in the supplied experiment scripts.
% RC_TLR32 denotes the full-grid spline input route.
    if nargin<1, options=struct; end
    c=struct('nValues',[256 512 1024 2048],'NValues',100:100:600, ...
        'pValues',[9 17 33 65 129 257 513 1025],'p',129, ...
        'qrTol',4e-8,'target',1e-7,'strength',2, ...
        'repeats',3,'warmups',1,'rankRepeats',1,'workspaceMB',128, ...
        'runLattice',true,'includeAlgorithm32',true,'runRankStudy',true,'runInterpolation',true, ...
        'plots',true,'plotPoints',161,'resume',true, ...
        'outDir',fullfile(fileparts(mfilename('fullpath')),'results'));
    names=fieldnames(options);
    for k=1:numel(names)
        assert(isfield(c,names{k}),'Unknown option: %s',names{k});
        c.(names{k})=options.(names{k});
    end
    jobs=struct('id',{},'system',{},'n',{},'N',{},'p',{},'method',{}, ...
        'interp',{},'filename',{},'protocol',{},'tol',{});
    for n=c.nValues(:).'
        add('biomol',n,500,c.p,'spline','RC_TLR');
    end
    if c.runRankStudy
        for N=c.NValues(:).'
            % Include p explicitly for the storage and coefficient spectra.
            for p=unique([c.pValues(:).' c.p],'stable')
                add('biomol',2048,N,p,'spline','RC_TLR');
            end
        end
    end
    if c.runInterpolation
        for type={'spline','linear','nearest'}
            add('biomol',2048,500,c.p,type{1},'RC_TLR');
        end
    end
    if c.includeAlgorithm32
        add('biomol',256,500,c.p,'spline','RC_TLR32');
    end
    if c.runLattice
        add('lattice',960,2304,c.p,'spline','RC_TLR');
    end
    [~,i]=unique({jobs.id},'stable'); jobs=jobs(i);
    [~,i]=sort(string({jobs.filename})); jobs=jobs(i);
    function add(system,n,N,p,interp,method)
        if strcmp(system,'biomol'), file=sprintf('n%d_Pn%d.mat',n,N);
        else, file='Lattice_n960_ParticlesNumber_2304_Rank_50.mat'; end
        id=sprintf('%s_n%d_N%d_p%d_%s_%s',system,n,N,p,method,interp);
        jobs(end+1)=struct('id',id,'system',system,'n',n,'N',N,'p',p, ...
            'method',method,'interp',interp,'filename',file, ...
            'protocol','standard','tol',c.qrTol); %#ok<AGROW>
    end
end
