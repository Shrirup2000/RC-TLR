function out = rc_tlr_biomol_experiment(dataRoot,options)
%RC_TLR_BIOMOL_EXPERIMENT Standalone RC-TLR biomolecular and optional lattice studies.
% out=rc_tlr_biomol_experiment('C:\project\data');
% Requires base MATLAB and the private helpers included in this folder.
% p is a NODE COUNT (degree <= p-1).
% Validation is the complete middle slice, not the full 3-D tensor.
% All times exclude input loading, CP references, diagnostics and plotting.
% TotalTime equals BuildSeconds. Biomolecular RC-TLR uses one timed build
% of the core, extension matrices and grid Lagrange bases. Other cases time
% continuous construction. Factor products and slice reconstruction are separate.
% Checkpoint resumes completed cases; failed/missing cases remain explicit.
    here=fileparts(mfilename('fullpath'));
    if nargin<1 || isempty(dataRoot), dataRoot=fullfile(here,'data'); end
    if nargin<2 || isempty(options), options=struct; end
    [jobs,cfg]=biomol_experiment_plan(options);
    validate_cfg(cfg); dataRoot=char(dataRoot);
    assert(isfolder(dataRoot),'Data folder does not exist: %s',dataRoot);
    if isfile(fullfile(dataRoot,'manifest.csv')) && ...
            any(~cellfun(@(f)isfile(fullfile(dataRoot,f)),unique({jobs.filename})))
        restore_biomol_data(dataRoot,unique({jobs.filename}));
    end
    if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
    files=dir(fullfile(dataRoot,'**','*.mat'));
    signature=struct('Version','rc-tlr-biomol-1.0','MATLAB',version, ...
        'DataRoot',dataRoot,'Settings',rmfield(cfg,{'outDir','resume'}));
    signature.TimeDefinition = ['Biomolecular RC-TLR: single factored grid build. ', ...
        'Other cases: continuous construction.'];
    signature.Environment=struct('Computer',computer,'BLAS',version('-blas'));
    % Changing any construction helper invalidates a resumed timing run.
    code=[dir(fullfile(here,'*.m'));dir(fullfile(here,'private','*.m'))];
    signature.Code=cell(numel(code),3);
    for k=1:numel(code)
        signature.Code(k,:)={fullfile(code(k).folder,code(k).name),code(k).bytes,code(k).datenum};
    end
    rows=repmat(blank_row(),0,1); raw=repmat(blank_time(),0,1);
    checkpoint=fullfile(cfg.outDir,'checkpoint.mat');
    if cfg.resume && isfile(checkpoint)
        old=load(checkpoint,'rows','raw','signature');
        assert(isequaln(old.signature,signature), ...
            'Settings/code changed. Use a new outDir or set resume=false.');
        rows=old.rows; raw=old.raw;
    end
    writetable(struct2table(jobs),fullfile(cfg.outDir,'experiment_plan.csv'));
    active=''; CP=[]; grids=[]; reference=[];
    fprintf('\nRC-TLR: %d configurations\n',numel(jobs));
    fprintf('Validation: complete middle slice. p counts nodes.\n');
    fprintf('Strong RRQR tolerance %.1e; final error target %.1e.\n',cfg.qrTol,cfg.target);
    for a=1:numel(jobs)
        j=jobs(a); row=blank_row(); row.CaseID=string(j.id);
        row.System=string(j.system); row.n=j.n; row.N=j.N; row.p=j.p;
        row.Method=string(j.method); row.Interpolation=string(j.interp);
        row.Protocol=string(j.protocol); row.Tolerance=j.tol;
        prev=[];
        if ~isempty(rows), prev=find([rows.CaseID]==row.CaseID,1); end
        try
            path=resolve_file(files,j.filename,dataRoot);
            s=dir(path); row.SourceFile=string(path);
            row.SourceBytes=s.bytes; row.SourceDatenum=s.datenum;
            if ~isempty(prev) && rows(prev).Status=="OK" && ...
                    rows(prev).SourceFile==row.SourceFile && ...
                    rows(prev).SourceBytes==s.bytes && rows(prev).SourceDatenum==s.datenum && ...
                    artifacts_present(j,cfg)
                fprintf('[%d/%d] Cached %s\n',a,numel(jobs),method_label(j.id)); continue
            end
            fprintf('[%d/%d] %s\n',a,numel(jobs),method_label(j.id));
            if ~strcmp(active,path)
                active=''; CP=[]; reference=[]; grids=[];
                [CP,grids]=load_input(path,j.system);
                validate_input(CP,grids,j.n);
                reference=cp_slice(CP,max(1,floor(j.n/2)));
                active=path;
            end
            dense=[];
            if strcmp(j.method,'RC_TLR32')
                assert(j.n<=256,'Algorithm 3.2 is restricted to n<=256 for 8 GB RAM.');
                assert(j.p<=129,'For the 8 GB driver, Algorithm 3.2 requires p<=129.');
                % Input preparation is outside the paper's approximation timer.
                dense=zeros(j.n,j.n,j.n);
                for k=1:j.n, dense(:,:,k)=cp_slice(CP,k); end
            end
            important=strcmp(j.system,'lattice') || (j.N==500 && j.p==cfg.p);
            if important, reps=cfg.repeats; warms=cfg.warmups;
            else, reps=cfg.rankRepeats; warms=0; end
            if is_biomol_rc_case(j)
                reps=1;
                warms=0;
            end
            times=zeros(reps,4); firstRanks=[]; firstSlice=[];
            if ~isempty(raw), raw([raw.CaseID]==row.CaseID)=[]; end
            for rep=1:(warms+reps)
                [M,build]=construct(CP,grids,j,cfg,dense);
                t=tic; V=evaluate_factors(M,grids); factor=toc(t);
                t=tic; Y=tucker_slice(M.C,V,max(1,floor(j.n/2))); recon=toc(t);
                rr=[size(M.C,1),size(M.C,2),size(M.C,3)];
                h=blank_time(); h.CaseID=row.CaseID; h.IsWarmup=rep<=warms;
                h.Repetition=max(0,rep-warms); h.BuildSeconds=build;
                h.FactorSeconds=factor; h.ReconstructionSeconds=recon;
                h.TotalSeconds=build; raw(end+1,1)=h; %#ok<AGROW>
                if rep<=warms, clear M V Y; continue; end
                times(rep-warms,:)=[build,factor,recon,build];
                if rep==warms+1
                    firstRanks=rr; firstSlice=Y;
                    row=measure_row(row,M,CP,Y,reference,grids,cfg);
                    if cfg.plots && important
                        save_surfaces(M,CP,Y,reference,grids,j,cfg);
                    end
                    if j.n==2048 && j.p==cfg.p && strcmp(j.interp,'spline') && ...
                            strcmp(j.method,'RC_TLR')
                        save_singular_values(M,j,cfg);
                    end
                else
                    assert(isequal(rr,firstRanks) && ...
                        norm(Y(:)-firstSlice(:))<=1e-10*max(norm(firstSlice(:)),realmin), ...
                        'Rebuilt approximation changed; timing repetitions cannot share one error.');
                end
                clear M V Y
            end
            if reps==1
                elapsed=times(1,:);
            else
                elapsed=median(times,1);
            end
            row.BuildSeconds=elapsed(1);
            row.FactorSeconds=elapsed(2);
            row.ReconstructionSeconds=elapsed(3);
            row.TotalTime=row.BuildSeconds;
            row.TotalTimeIQR=quartile(times(:,1),.75)-quartile(times(:,1),.25);
            row.Repeats=reps; row.Status="OK";
            fprintf('  ranks [%d %d %d], abs %.3e, rel %.3e, total %.4g s\n', ...
                row.rX,row.rY,row.rZ,row.AbsoluteError,row.RelativeError,row.TotalTime);
            clear dense firstSlice
        catch ME
            row.Status="FAILED"; row.Message=string(ME.message);
            row.ErrorIdentifier=string(ME.identifier);
            report=getReport(ME,'extended','hyperlinks','off');
            fid=fopen(fullfile(cfg.outDir,[j.id '_FAILED.txt']),'w');
            if fid~=-1, fprintf(fid,'%s\n',report); fclose(fid); end
            fprintf('  FAILED: %s\n',ME.message);
            clear M V Y dense firstSlice
        end
        if isempty(prev), rows(end+1,1)=row; else, rows(prev)=row; end %#ok<AGROW>
        tmp=[tempname(cfg.outDir) '.mat'];
        save(tmp,'rows','raw','signature','-v7'); movefile(tmp,checkpoint,'f');
        writetable(report_table(struct2table(rows)),fullfile(cfg.outDir,'all_cases.csv'));
    end
    out=struct('Settings',cfg,'Metadata',signature,'Cases',struct2table(rows), ...
        'Repetitions',struct2table(raw),'Plan',struct2table(jobs));
    out.AllRequestedCasesOK=height(out.Cases)==numel(jobs) && all(out.Cases.Status=="OK");
    out.Tables=make_tables(out.Cases,cfg);
    names=fieldnames(out.Tables);
    for k=1:numel(names)
        writetable(report_table(out.Tables.(names{k})),fullfile(cfg.outDir,[names{k} '.csv']));
    end
    writetable(out.Repetitions,fullfile(cfg.outDir,'timing_repetitions.csv'));
    out.PlotError="";
    if cfg.plots
        try, summary_plots(out.Cases,cfg);
        catch ME
            out.PlotError=string(getReport(ME,'extended','hyperlinks','off'));
            warning('Summary plot: %s',ME.message);
        end
    end
    save(fullfile(cfg.outDir,'results.mat'),'out','-v7');
    fprintf('\nCompleted %d / %d cases. Results: %s\n', ...
        nnz(out.Cases.Status=="OK"),height(out.Cases),cfg.outDir);
    if any(out.Cases.Status~="OK")
        disp(out.Cases(out.Cases.Status~="OK",{'CaseID','Message'}));
    end
    show={'n','N','Method','Tolerance','CoreStorage','AbsoluteError','RelativeError', ...
        'ContinuousStorage','GridStorage','FactoredGridStorage','BuildSeconds','TotalTime','Status'};
    fprintf('\nGrid-size study\n'); disp(report_table(out.Tables.Table3(:,show)));
    fprintf('\nTable 4: RC-TLR maximum ranks\n'); disp(out.Tables.Ranks_RC_TLR);
    fprintf('\nStorage study\n');
    showTable5={'N','Method','CoreStorage','CPFactorStorage', ...
        'ContinuousStorage','CompressionPct','SavingsPct','Status'};
    disp(report_table(out.Tables.Table5(:,showTable5)));
    if cfg.runLattice, fprintf('\nLattice study\n'); disp(report_table(out.Tables.Lattice(:,show))); end
end

function validate_cfg(c)
    for name={'nValues','NValues','pValues'}
        validateattributes(c.(name{1}),{'numeric'},{'vector','integer','positive','finite','nonempty'});
    end
    for name={'p','repeats','rankRepeats','plotPoints'}
        validateattributes(c.(name{1}),{'numeric'},{'scalar','integer','positive','finite'});
    end
    validateattributes(c.warmups,{'numeric'},{'scalar','integer','nonnegative','finite'});
    validateattributes(c.workspaceMB,{'numeric'},{'scalar','>=',16,'finite'});
    for name={'qrTol','target'}
        validateattributes(c.(name{1}),{'numeric'},{'scalar','>',0,'<',1,'finite'});
    end
    validateattributes(c.strength,{'numeric'},{'scalar','>',1,'finite'});
end

function path=resolve_file(files,name,root)
    k=find(strcmp({files.name},name));
    assert(~isempty(k),'Missing input %s under %s.',name,root);
    assert(numel(k)==1,'Multiple copies of %s. Point dataRoot at one dataset tree.',name);
    path=fullfile(files(k).folder,files(k).name);
end

function [CP,g]=load_input(path,system)
    if strcmp(system,'biomol')
        S=load(path,'LRed','Hunif');
        CP=struct('lambda',S.LRed.LAM1C(:)./S.Hunif^3, ...
            'U',{{S.LRed.CU1,S.LRed.CU2,S.LRed.CU3}});
        CP.originalRank=numel(CP.lambda); CP.removedShortRank=0;
        CP.vectorScale=1/S.Hunif;
        g=cellfun(@(u)linspace(-1+1/size(u,1),1-1/size(u,1),size(u,1)), ...
            CP.U,'UniformOutput',false);
    else
        S=load(path,'Ldip','Lat','ASC','AvC','xcol','ycol','zcol');
        Rp=size(S.ASC.U1k,2); Rn=size(S.AvC.U1k,2); Rt=Rp+Rn;
        keep=[1:Rp-11,Rp+1:Rt-11]; % Exactly the supplied lattice_plots split.
        U={S.Ldip.Ut1./(2/S.Lat.x(end)),S.Ldip.Ut2./(2/S.Lat.y(end)), ...
            S.Ldip.Ut3./(1/S.Lat.z(end))};
        w=[ones(Rp,1);-ones(Rn,1)];
        U=cellfun(@(a)a(:,keep),U,'UniformOutput',false);
        Q=CP_normalize({w(keep),U});
        CP=struct('lambda',Q{1}(:),'U',{Q{2}});
        CP.originalRank=Rt; CP.removedShortRank=Rt-numel(keep);
        CP.vectorScale=1;
        g={S.xcol(:).',S.ycol(:).',S.zcol(:).'};
    end
    CP.U=reshape(CP.U,1,3); g=reshape(g,1,3);
end

function validate_input(CP,g,n)
    R=numel(CP.lambda); assert(all(isfinite(CP.lambda)),'Nonfinite CP weights.');
    for d=1:3
        validateattributes(CP.U{d},{'numeric'},{'size',[n R],'real','finite'});
        expected=linspace(-1+1/n,1-1/n,n);
        assert(numel(g{d})==n && max(abs(g{d}(:).'-expected))<1e-11, ...
            'This experiment uses cell-centred coordinates on [-1,1].');
    end
end

function [M,seconds]=construct(CP,g,j,c,dense)
    p=j.p; domain=repmat([-1 1],3,1); t=tic;
    M=struct('method',j.method,'kind','','C',[],'H',{{}},'raw',[], ...
        'nodes',my_chev2(p),'weights',CP.lambda, ...
        'tol',NaN,'qrResidual',NaN,'rho',NaN);
    if strcmp(j.method,'RC_TLR')
        if strcmp(j.interp,'spline')
            [Gw,~,~,~,G]=weighted_cheb_tensor(p,CP,g,domain);
        else
            P=cell(1,3);
            for d=1:3
                P{d}=interp1(g{d}(:),CP.U{d},M.nodes(:),j.interp,'extrap');
            end
            G=struct('lambda',CP.lambda,'U',{P},'weighted',false);
            [wx,wy,wz]=cheb_weight_matrices(p); Gw=weighted_core(G,wx,wy,wz);
        end
        if is_biomol_rc_case(j) && p == 129
            % Assemble the weighted samples once for the three mode unfoldings.
            GwDense = rc_tlr_cp_values(Gw);
            [I,~,~,q] = select_cheb_nodes(GwDense,j.tol,c.strength);
            for d = 1:3
                q{d}.weights = Gw.weights{d};
                assert(q{d}.accuracySatisfied && q{d}.strongSatisfied, ...
                    'RC-TLR node selection failed its QR checks in mode %d.',d);
            end
            clear GwDense
        else
            [I,~,q]=rc_tlr_bounded_nodes(Gw,j.tol,c.strength,c.workspaceMB);
        end
        M.C=selected_core(p,I,G);
        [~,~,~,M.H]=selected_factor_matrices({[],[],[]},p,I,q,domain);
        M.kind='values'; M.raw=G.U; M.tol=j.tol;
        M.qrResidual=max(cellfun(@(v)v.relativeResidual,q));
        M.rho=max(cellfun(@(v)v.maxRho,q));
    elseif strcmp(j.method,'RC_TLR32')
        % Algorithm 2: full tensor -> multivariate spline -> functional RC.
        opts=struct('strength',c.strength,'maxSampleMB',128);
        [F,G]=rc_tlr_from_grid(dense,g,[p p p],j.tol,opts);
        M.C=F.C; M.H=F.H; M.raw=G; M.nodes=F.nodes{1};
        M.kind='dense_values'; M.tol=j.tol;
        M.qrResidual=F.MaxQRResidual; M.rho=F.MaxRho;
    else
        error('Unsupported RC-TLR input route: %s',j.method);
    end
    if is_biomol_rc_case(j)
        M.B=cell(1,3);
        for d=1:3
            M.B{d}=lagrange_cheb(M.nodes,g{d});
        end
    end
    seconds=toc(t);
end

function yes=is_biomol_rc_case(j)
    yes=strcmp(j.system,'biomol') && strcmp(j.method,'RC_TLR');
end

function V=evaluate_factors(M,g)
    V=cell(1,3);
    for d=1:3
        if isfield(M,'B')
            V{d}=M.B{d}*M.H{d};
        else
            V{d}=bary_eval(M.nodes,M.H{d},g{d});
        end
    end
end

function V=bary_eval(nodes,H,x)
    V=zeros(numel(x),size(H,2));
    for first=1:128:numel(x)
        ix=first:min(first+127,numel(x));
        V(ix,:)=lagrange_cheb(nodes,x(ix))*H;
    end
end

function Y=tucker_slice(C,V,k)
    r=[size(C,1),size(C,2),size(C,3)];
    B=reshape(reshape(C,r(1)*r(2),r(3))*V{3}(k,:).',r(1),r(2));
    Y=(V{1}*B)*V{2}.';
end

function Y=cp_slice(CP,k)
    Y=(CP.U{1}.*(CP.lambda(:).'.*CP.U{3}(k,:)))*CP.U{2}.';
end

function Y=uncompressed_slice(M,g)
    k=max(1,floor(numel(g{3})/2));
    if strcmp(M.kind,'dense_values')
        V={lagrange_cheb(M.nodes,g{1}),lagrange_cheb(M.nodes,g{2}), ...
            lagrange_cheb(M.nodes,g{3}(k))};
        Y=tucker_slice(M.raw,V,1);
    else
        V=cell(1,3);
        for d=1:3
            if d==3, x=g{d}(k); else, x=g{d}; end
            V{d}=bary_eval(M.nodes,M.raw{d},x);
        end
        Y=(V{1}.*(M.weights(:).'.*V{3}))*V{2}.';
    end
end

function row=measure_row(row,M,CP,Y,reference,g,c)
    rr=[size(M.C,1),size(M.C,2),size(M.C,3)];
    row.Scope="middle_slice"; row.MiddleIndex=max(1,floor(row.n/2));
    row.ValidationEntries=numel(Y);
    row.PolynomialDegree=row.p-1;
    row.NodeFamily="first_kind"; row.CompressionBackend="strong_RRQR";
    row.InputRank=numel(CP.lambda); row.OriginalCPRank=CP.originalRank;
    row.SourceCPRank=numel(CP.lambda); row.InputRepresentation="CP";
    if row.Method=="RC_TLR32"
        row.InputRepresentation="full_grid"; row.InputRank=NaN;
    end
    row.RemovedShortRank=CP.removedShortRank;
    row.rX=rr(1); row.rY=rr(2); row.rZ=rr(3);
    row.CoreStorage=string(sprintf('%d x %d x %d',rr));
    row.MaxRank=max(rr); row.CoreEntries=numel(M.C);
    row.ContinuousStorage=numel(M.C)+sum(cellfun(@numel,M.H));
    if isfield(M,'B')
        row.FactoredGridStorage=row.ContinuousStorage+sum(cellfun(@numel,M.B));
    end
    row.GridStorage=numel(M.C)+sum(cellfun(@numel,g).*rr);
    row.CPFactorStorage=sum(cellfun(@numel,CP.U));
    row.CPStorage=row.CPFactorStorage+numel(CP.lambda);
    row.InputStorage=row.CPStorage;
    if row.InputRepresentation=="full_grid", row.InputStorage=row.n^3; end
    row.InputCompressionPct=100*row.ContinuousStorage/row.InputStorage;
    row.InputSavingsPct=100-row.InputCompressionPct;
    row.ContinuousMB=8*row.ContinuousStorage/2^20;
    row.GridMB=8*row.GridStorage/2^20;
    row.CompressionPct=100*row.ContinuousStorage/row.CPFactorStorage;
    row.SavingsPct=100-row.CompressionPct;
    delta=Y-reference; row.AbsoluteError=max(abs(delta(:)));
    row.RelativeError=ratio(row.AbsoluteError,max(abs(reference(:))));
    row.RelativeFrobeniusError=ratio(norm(delta(:)),norm(reference(:)));
    P=uncompressed_slice(M,g);
    row.CompressionError=max(abs(P(:)-Y(:)));
    row.InterpolationError=max(abs(P(:)-reference(:)));
    row.AbsoluteToleranceMet=row.AbsoluteError<=c.target;
    row.RelativeToleranceMet=row.RelativeError<=c.target;
    row.Tolerance=M.tol; row.MaxQRResidual=M.qrResidual; row.MaxRho=M.rho;
end

function save_surfaces(M,CP,Y,reference,g,j,c) %#ok<INUSD>
    n=numel(g{1}); ix=unique(round(linspace(1,n,min(n,c.plotPoints))));
    P=uncompressed_slice(M,g);
    % Signed errors follow approximate minus reference conventions.
    data={Y(ix,ix),Y(ix,ix)-P(ix,ix),Y(ix,ix)-reference(ix,ix)};
    titles={'Approximation','Compression error','Total error'};
    f=figure('Visible','off','Color','w','Position',[100 100 1200 350]);
    cleanup=onCleanup(@()close(f)); %#ok<NASGU>
    for k=1:3
        subplot(1,3,k); surf(g{1}(ix),g{2}(ix),data{k}.','EdgeColor','none');
        title(titles{k}); xlabel('x'); ylabel('y'); axis tight; view(3); colorbar;
    end
    sgtitle(method_label(j.id));
    exportgraphics(f,fullfile(c.outDir,[j.id '.png']),'Resolution',130);
    % Retain only plotting samples, not complete slices or CP inputs.
    x=g{1}(ix); y=g{2}(ix); Ref=reference(ix,ix); Approx=data{1}; %#ok<NASGU>
    Compression=data{2}; Total=data{3}; %#ok<NASGU>
    save(fullfile(c.outDir,[j.id '_slice.mat']),'x','y','Ref','Approx','Compression','Total');
    % Original long-range input and ten x-direction CP factors (paper Fig. 1).
    if strcmp(j.method,'RC_TLR') && strcmp(j.interp,'spline')
        clear cleanup
        f=figure('Visible','off','Color','w','Position',[100 100 900 350]);
        cleanup=onCleanup(@()close(f)); %#ok<NASGU>
        subplot(1,2,1); surf(x,y,Ref.','EdgeColor','none');
        title('Input long-range CP slice'); xlabel('x'); ylabel('y'); view(3); axis tight
        factors=CP.U{1}(:,1:min(10,numel(CP.lambda)))*CP.vectorScale;
        subplot(1,2,2); plot(g{1},factors);
        title('First ten x-direction CP factors'); xlabel('x'); grid on
        exportgraphics(f,fullfile(c.outDir,[j.id '_input.png']),'Resolution',130);
        if strcmp(j.system,'biomol') && j.n==1024
            clear cleanup
            f=figure('Visible','off','Color','w'); cleanup=onCleanup(@()close(f));
            plot(g{1},factors); xlabel('x'); ylabel('Canonical factor'); grid on
            title('Figure 1: first ten biomolecular x factors, n = 1024');
            exportgraphics(f,fullfile(c.outDir,'Figure1_biomol_vectors.png'),'Resolution',150);
            writematrix([g{1}(:) factors],fullfile(c.outDir,'Figure1_biomol_vectors.csv'));
            clear cleanup
            f=figure('Visible','off','Color','w'); cleanup=onCleanup(@()close(f));
            surf(x,y,Ref.','EdgeColor','none'); xlabel('x'); ylabel('y'); view(3); axis tight
            title(sprintf('Figure 4 left: supplied long-range input, rank %d',numel(CP.lambda)));
            exportgraphics(f,fullfile(c.outDir,'Figure4_left_input.png'),'Resolution',150);
        end
    end
end

function save_singular_values(M,j,c)
    K=cell(1,3);
    for d=1:3
        if strcmp(M.kind,'values')
            p=j.p; theta=(2*(1:p)-1)*pi/(2*p);
            K{d}=(2/p)*cos((0:p-1).'*theta)*M.raw{d}; K{d}(1,:)=K{d}(1,:)/2;
        else, K{d}=M.raw{d}; end
    end
    v=cellfun(@(a)svd(a,'econ'),K,'UniformOutput',false);
    S=table((1:numel(v{1})).',v{1},v{2},v{3}, ...
        'VariableNames',{'Index','Mode1','Mode2','Mode3'});
    writetable(S,fullfile(c.outDir,[j.id '_singular_values.csv']));
end

function T=make_tables(C,c)
    bio=C.System=="biomol"; spline=C.Interpolation=="spline";
    grid=bio & C.N==500 & C.p==c.p & spline & ismember(C.n,c.nValues) & C.Method=="RC_TLR";
    standard=C.Protocol=="standard";
    T.Table3=sortrows(C(grid & standard,:),{'n','Method'});
    T.Algorithm32=C(bio & C.Method=="RC_TLR32",:);
    T.Table4=sortrows(C(bio & C.n==2048 & spline & standard & C.Method=="RC_TLR" & ...
        ismember(C.N,c.NValues) & ismember(C.p,c.pValues),:),{'N','p','Method'});
    T.Table5=sortrows(C(bio & C.n==2048 & spline & standard & C.p==c.p & ...
        ismember(C.N,c.NValues) & C.Method=="RC_TLR",:),{'N','Method'});
    T.Interpolation=C(bio & C.n==2048 & C.N==500 & C.p==c.p & standard & C.Method=="RC_TLR",:);
    T.Lattice=C(C.System=="lattice",:);
    A=table(c.NValues(:),'VariableNames',{'N'});
    for p=c.pValues(:).'
        v=nan(numel(c.NValues),1);
        for k=1:numel(v)
            r=T.Table4(T.Table4.N==c.NValues(k) & ...
                T.Table4.p==p & T.Table4.Status=="OK",:);
            if ~isempty(r), v(k)=r.MaxRank(1); end
        end
        A.(sprintf('p%d',p))=v;
    end
    T.Ranks_RC_TLR=A;
end

function summary_plots(C,c)
    C=C(C.Status=="OK" & C.Protocol=="standard",:);
    allCases=C; C=C(C.Interpolation=="spline",:);
    fields={'AbsoluteError','RelativeError','ContinuousStorage','TotalTime'};
    labels={'Middle-slice absolute maximum error','Middle-slice relative maximum error', ...
        'Continuous parameters','Construction time (s)'};
    method="RC_TLR";
    f=figure('Visible','off','Color','w','Position',[100 100 1100 750]);
    cleanup=onCleanup(@()close(f)); %#ok<NASGU>
    for k=1:4
        subplot(2,2,k); hold on
        A=C(C.System=="biomol" & C.N==500 & C.p==c.p & C.Method==method & ismember(C.n,c.nValues),:);
        A=sortrows(A,'n');
        semilogy(A.n,max(A.(fields{k}),realmin),'-o','DisplayName',method_label(method));
        set(gca,'YScale','log'); xlabel('Grid points per coordinate n'); ylabel(labels{k});
        grid on; legend('Location','best');
    end
    exportgraphics(f,fullfile(c.outDir,'grid_study.png'),'Resolution',150);
    clear cleanup
    f=figure('Visible','off','Color','w','Position',[100 100 1000 400]);
    cleanup=onCleanup(@()close(f)); %#ok<NASGU>
    for k=1:2
        subplot(1,2,k); hold on
        A=C(C.System=="biomol" & C.n==2048 & C.p==c.p & C.Method==method & ismember(C.N,c.NValues),:);
        A=sortrows(A,'N');
        if k==1, plot(A.N,A.MaxRank,'-o','DisplayName',method_label(method));
        else, plot(A.N,A.ContinuousStorage,'-o','DisplayName',method_label(method)); end
        xlabel('Atoms N'); grid on; legend('Location','best');
        if k==1, ylabel('Maximum Tucker rank'); else, ylabel('Continuous parameters'); end
    end
    exportgraphics(f,fullfile(c.outDir,'rank_storage_scaling.png'),'Resolution',150);
    clear cleanup
    % Figure 7: interpolation methods on a common complete reference slice.
    for method="RC_TLR"
        f=figure('Visible','off','Color','w','Position',[100 100 1200 350]);
        cleanup=onCleanup(@()close(f)); %#ok<NASGU>
        types={'spline','linear','nearest'};
        for k=1:3
            name=sprintf('biomol_n2048_N500_p%d_%s_%s_slice.mat',c.p,method,types{k});
            path=fullfile(c.outDir,name);
            if ~isfile(path), continue; end
            S=load(path); subplot(1,3,k);
            surf(S.x,S.y,S.Total.','EdgeColor','none');
            title(types{k}); xlabel('x'); ylabel('y'); axis tight; view(3); colorbar
        end
        sgtitle(method_label(method));
        exportgraphics(f,fullfile(c.outDir,[char(method) '_interpolation.png']),'Resolution',130);
        if c.p==129 && has_cases(allCases,method,[2048 2048 2048],500,129, ...
                ["spline" "linear" "nearest"])
            copyfile(fullfile(c.outDir,[char(method) '_interpolation.png']), ...
                fullfile(c.outDir,['Figure7_' char(method) '.png']));
        end
        clear cleanup
        % All three side-matrix spectra are exported; plot mode 2 as Fig. 8.
        f=figure('Visible','off','Color','w','Position',[100 100 1000 400]);
        cleanup=onCleanup(@()close(f)); %#ok<NASGU>
        for N=c.NValues(:).'
            name=sprintf('biomol_n2048_N%d_p%d_%s_spline_singular_values.csv',N,c.p,method);
            path=fullfile(c.outDir,name);
            if ~isfile(path), continue; end
            S=readtable(path);
            for k=1:2
                subplot(1,2,k); hold on
                semilogy(S.Index,max(S.Mode2,realmin),'DisplayName',sprintf('N = %d',N));
                set(gca,'YScale','log'); xlabel('Singular-value index'); ylabel('Singular value');
                if k==2, xlim([10 40]); end
                grid on; legend('Location','best');
            end
        end
        sgtitle([method_label(char(method)) ': mode-2 coefficient-factor spectra']);
        exportgraphics(f,fullfile(c.outDir,[char(method) '_singular_values.png']),'Resolution',150);
        if c.p==129 && has_cases(C,method,2048*ones(1,6),100:100:600,129,"spline")
            copyfile(fullfile(c.outDir,[char(method) '_singular_values.png']), ...
                fullfile(c.outDir,['Figure8_' char(method) '.png']));
        end
        clear cleanup
        if c.p==129
            % Figures 5/9 use the already exported approximation/error triplets.
            for spec={ {5,'biomol',256,500}, {9,'lattice',960,2304} }
                s=spec{1}; name=sprintf('%s_n%d_N%d_p129_%s_spline.png',s{2},s{3},s{4},method);
                if has_cases(C,method,s{3},s{4},129,"spline") && isfile(fullfile(c.outDir,name))
                    copyfile(fullfile(c.outDir,name), ...
                        fullfile(c.outDir,sprintf('Figure%d_%s.png',s{1},method)));
                end
            end
            if has_cases(C,method,[512 1024 2048],500,129,"spline")
                f=figure('Visible','off','Color','w','Position',[100 100 1200 350]);
                cleanup=onCleanup(@()close(f)); %#ok<NASGU>
                ns=[512 1024 2048];
                for k=1:3
                    name=sprintf('biomol_n%d_N500_p129_%s_spline_slice.mat',ns(k),method);
                    S=load(fullfile(c.outDir,name)); subplot(1,3,k);
                    surf(S.x,S.y,S.Total.','EdgeColor','none');
                    title(sprintf('n = %d',ns(k))); xlabel('x'); ylabel('y'); view(3); axis tight
                end
                sgtitle(['Figure 6: total errors, ' method_label(char(method))]);
                exportgraphics(f,fullfile(c.outDir,['Figure6_' char(method) '.png']),'Resolution',150);
                clear cleanup
            end
        end
    end
end

function yes=has_cases(C,method,ns,Ns,p,types)
    count=max([numel(ns),numel(Ns),numel(types)]); yes=true;
    for i=1:count
        n=ns(min(i,numel(ns))); N=Ns(min(i,numel(Ns))); type=types(min(i,numel(types)));
        yes=yes && any(C.Method==method & C.n==n & C.N==N & C.p==p & C.Interpolation==type);
    end
end

function yes=artifacts_present(j,c)
    yes=true;
    important=strcmp(j.system,'lattice') || (j.N==500 && j.p==c.p);
    if c.plots && important
        yes=isfile(fullfile(c.outDir,[j.id '.png'])) && ...
            isfile(fullfile(c.outDir,[j.id '_slice.mat']));
        if strcmp(j.method,'RC_TLR') && strcmp(j.system,'biomol') && j.n==1024
            yes=yes && isfile(fullfile(c.outDir,'Figure1_biomol_vectors.png')) && ...
                isfile(fullfile(c.outDir,'Figure4_left_input.png'));
        end
    end
    if j.n==2048 && j.p==c.p && strcmp(j.interp,'spline') && strcmp(j.method,'RC_TLR')
        yes=yes && isfile(fullfile(c.outDir,[j.id '_singular_values.csv']));
    end
end

function v=ratio(a,b)
    if b==0, if a==0, v=0; else, v=Inf; end
    else, v=a/b; end
end

function v=quartile(x,q)
    x=sort(x(:)); h=numel(x)*q+.5;
    if h<=1, v=x(1); elseif h>=numel(x), v=x(end);
    else, k=floor(h); v=x(k)+(h-k)*(x(k+1)-x(k)); end
end

function r=blank_row()
    r=struct('CaseID',"",'System',"",'Method',"",'Interpolation',"",'Protocol',"", ...
        'Scope',"middle_slice",'MiddleIndex',NaN,'ValidationEntries',NaN, ...
        'n',NaN,'N',NaN,'p',NaN,'PolynomialDegree',NaN,'NodeFamily',"",'CompressionBackend',"", ...
        'InputRank',NaN,'SourceCPRank',NaN,'InputRepresentation',"", ...
        'OriginalCPRank',NaN,'RemovedShortRank',NaN,'InputStorage',NaN, ...
        'InputCompressionPct',NaN,'InputSavingsPct',NaN, ...
        'rX',NaN,'rY',NaN,'rZ',NaN,'CoreStorage',"", ...
        'MaxRank',NaN,'CoreEntries',NaN,'ContinuousStorage',NaN,'GridStorage',NaN, ...
        'FactoredGridStorage',NaN, ...
        'CPFactorStorage',NaN,'CPStorage',NaN,'CompressionPct',NaN,'SavingsPct',NaN, ...
        'ContinuousMB',NaN,'GridMB',NaN, ...
        'BuildSeconds',NaN,'FactorSeconds',NaN,'ReconstructionSeconds',NaN, ...
        'TotalTime',NaN,'TotalTimeIQR',NaN,'Repeats',NaN,'AbsoluteError',NaN, ...
        'RelativeError',NaN,'RelativeFrobeniusError',NaN,'CompressionError',NaN, ...
        'InterpolationError',NaN,'AbsoluteToleranceMet',false,'RelativeToleranceMet',false, ...
        'Tolerance',NaN,'MaxQRResidual',NaN,'MaxRho',NaN,'Status',"PENDING", ...
        'Message',"",'ErrorIdentifier',"",'SourceFile',"",'SourceBytes',NaN,'SourceDatenum',NaN);
end

function r=blank_time()
    r=struct('CaseID',"",'Repetition',0,'IsWarmup',false,'BuildSeconds',NaN, ...
        'FactorSeconds',NaN,'ReconstructionSeconds',NaN,'TotalSeconds',NaN);
end

function T=report_table(T)
    if any(strcmp('Method',T.Properties.VariableNames))
        T.Method=method_label(T.Method);
    end
end

function label=method_label(value)
    label=strrep(value,'RC_TLR32','RC-TLR (Alg. 3.2)');
    label=strrep(label,'RC_TLR','RC-TLR');
    label=strrep(label,'_',' ');
end
